/**
 * echidna.c
 *  - On-the-fly stream parallelizer for FASTQ and FASTA processing
 *
 * Copyright (c) 2012 Hyeshik Chang
 * 
 * Permission is hereby granted, free of charge, to any person obtaining a copy
 * of this software and associated documentation files (the "Software"), to deal
 * in the Software without restriction, including without limitation the rights
 * to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
 * copies of the Software, and to permit persons to whom the Software is
 * furnished to do so, subject to the following conditions:
 * 
 * The above copyright notice and this permission notice shall be included in
 * all copies or substantial portions of the Software.
 * 
 * THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
 * IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
 * FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE 
 * AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
 * LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
 * OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN
 * THE SOFTWARE.
 */

#define _GNU_SOURCE
#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <string.h>
#include <fcntl.h>
#include <signal.h>
#include <getopt.h>
#include <stdarg.h>
#include <err.h>
#include <errno.h>
#include <inttypes.h>
#include <sys/select.h>
#include <sys/types.h>
#include <sys/wait.h>
#include "bufqueue.h"
#include "config.h"

#define DEFAULT_NUM_WORKERS 4
#define INQUEUE_SIZE        65536
#define OUTQUEUE_SIZE       65536
#define INLET_QUEUE_SIZE    262144
#define OUTLET_QUEUE_SIZE   262144

#define STATUS_IDLE         0
#define STATUS_RUNNING      3
#define STATUS_FLUSHING     2
#define STATUS_TERMINATED   4

#define IS_STDIN_ALIVE(w)   ((w)->status & 1)
#define IS_STDOUT_ALIVE(w)  ((w)->status & 2)

#define FORMAT_UNDECIDED    0
#define FORMAT_FASTQ        1
#define FORMAT_FASTA        2

struct _WORKER;
struct _SESSION;

typedef struct _WORKER {
    int workerid;
    int status;
    pid_t pid;

    int stdin_fd;
    int stdout_fd;
    QUEUE *inbuf;
    QUEUE *outbuf;

    int (*input_handler)(struct _SESSION *sess, struct _WORKER *worker);
    uint64_t lineno;
} WORKER;

typedef struct _SESSION {
    int num_workers;
    int running_workers;
    WORKER *workers;

    int rr_next;

    QUEUE *inbuf;
    QUEUE *outbuf;

    int (*input_handler)(struct _SESSION *sess);
    uint64_t lineno;

    /* Only one between `command' and `args' which is non-NULL is effective */
    const char *command;
    char **args;
} SESSION;

static SESSION *global_session=NULL; /* for signal handlers */

#define FD_MAX(a, b) ((a) > (b) ? (a) : (b))

static void
error(int status, const char *msg, ...)
{
    va_list args;

    fprintf(stderr, "echidna: ");

    va_start(args, msg);
    vfprintf(stderr, msg, args);
    va_end(args);

    if (status != 0)
        exit(status);
}

static void
sigchld_hdl(int sig)
{
    int i;
    pid_t pid;

    while ((pid = waitpid(-1, NULL, WNOHANG)) > 0) {
        if (global_session->workers == NULL)
            continue;

        for (i = 0; i < global_session->num_workers; i++)
            if (global_session->workers[i].pid == pid) {
                global_session->workers[i].status = STATUS_TERMINATED;
                global_session->running_workers--;
                break;
            }
    }
}

static int
set_io_nonblocking(int fd)
{
    int flags = fcntl(fd, F_GETFL, 0);

    if (flags == -1)
        return -1;

    return fcntl(fd, F_SETFL, flags | O_NONBLOCK);
}

static int
handle_input_from_stdin_fastq(SESSION *sess)
{
    char *head, *tail, *cur, *corner;
    int i, line_in_record;

    head = cur = sess->inbuf->data + sess->inbuf->front;
    corner = sess->inbuf->data + sess->inbuf->size;
    tail = sess->inbuf->data + sess->inbuf->rear;

    line_in_record = 0;

    while (cur != tail) {
        if (*cur == '\n' && ++line_in_record == 4) {
            qsize_t recordsize;

            recordsize = (head <= cur ? cur + 1 - head :
                          (corner - head) + (cur + 1 - sess->inbuf->data));

            for (i = 0; i < sess->num_workers; i++) {
                int wselected=(sess->rr_next + i) % sess->num_workers;

                if (*head != '@')
                    error(0, "Unaligned FASTQ input at line %d\n",
                          sess->lineno);

                if (queue_transfer(sess->workers[wselected].outbuf,
                                   sess->inbuf, recordsize) != -1) {
                    sess->rr_next = (sess->rr_next+1) % sess->num_workers;
                    break;
                }
                /* TODO: handle records way too big in outbuf. */
                //errx(13, "line %ld: found a too large record. (1)\n",
                //         sess->lineno);
            }

            if (i == sess->num_workers) /* all buffers are full */
                break;

            // intended duplication of code for optimized running
            cur++;
            if (cur == corner)
                cur = sess->inbuf->data;

            head = cur;
            sess->lineno += 4;
            line_in_record = 0;

            if (head != sess->inbuf->data + sess->inbuf->front)
                error(0, "Internal programming error: head and front "
                         "mismatching - %p, %p, %d\n",
                         head, sess->inbuf->data, sess->inbuf->front);
        }
        else {
            cur++;
            if (cur == corner)
                cur = sess->inbuf->data;
        }
    }

    return 0;
}

static int
handle_input_from_stdin_undecided(SESSION *sess)
{
    if (queue_num_filled(sess->inbuf) >= 1)
        switch (sess->inbuf->data[sess->inbuf->front]) {
        case '@': /* FASTQ */
            sess->input_handler = handle_input_from_stdin_fastq;
            return handle_input_from_stdin_fastq(sess);
            break;
        case '>': /* FASTA */
            error(1, "FASTA support is not implemented yet.\n");
            break;
        default:
            error(1, "Unknown input format: the first character is not "
                     "'@' or '>'.\n");
        }

    return 0;
}

static int
handle_input_from_worker_fastq(SESSION *sess, WORKER *worker)
{
    char *head, *tail, *cur, *leftend, *rightend;
    int line_in_record;

    head = cur = worker->inbuf->data + worker->inbuf->front;
    leftend = worker->inbuf->data;
    rightend = worker->inbuf->data + worker->inbuf->size;
    tail = worker->inbuf->data + worker->inbuf->rear;

    line_in_record = 0;

    while (cur != tail) {
        if (*cur == '\n' && ++line_in_record == 4) {
            qsize_t recordsize;

            recordsize = (head <= cur ? cur + 1 - head :
                          (rightend - head) + (cur + 1 - leftend));

            if (*head != '@')
                error(0, "(worker %d) Unaligned FASTQ input at"
                        " line no %d\n", worker->workerid, worker->lineno);

            if (queue_transfer(sess->outbuf, worker->inbuf, recordsize) == -1)
                break;
            /* TODO: handle records way too big in outbuf. */
            //errx(13, "line %ld: found a too large record. (1)\n",
            //         sess->lineno);

            // intended duplication of code for optimized running
            cur++;
            if (cur == rightend)
                cur = leftend;

            head = cur;
            worker->lineno += 4;
            line_in_record = 0;

            if (head != leftend + worker->inbuf->front)
                error(0, "Internal programming error: head and front "
                         "mismatching - %p, %p, %d\n",
                         head, leftend, worker->inbuf->front);
        }
        else {
            cur++;
            if (cur == rightend)
                cur = leftend;
        }
    }

    return 0;
}

static int
handle_input_from_worker_fasta(SESSION *sess, WORKER *worker)
{
    char *head, *tail, *cur, *leftend, *rightend;
    qflag_t header_read;

    head = cur = worker->inbuf->data + worker->inbuf->front;
    leftend = worker->inbuf->data;
    rightend = worker->inbuf->data + worker->inbuf->size;
    tail = worker->inbuf->data + worker->inbuf->rear;
    header_read = worker->inbuf->flags;

    while (cur != tail) {
        if (*cur == '\n')
            worker->lineno++;

        if (!header_read) {
            if (*cur == '\n')
                header_read = 1;
        }
        else if (*cur == '>') {
            qsize_t recordsize;

            recordsize = (head <= cur ? cur - head : (rightend - head) + (cur - leftend));

            if (queue_transfer(sess->outbuf, worker->inbuf, recordsize) == -1)
                break; /* TODO: handle records way too big in outbuf. */

            head = cur;
            header_read = 0;
        }

        if (++cur == rightend)
            cur = leftend;
    }

    worker->inbuf->flags = header_read;
    return 0;
}

static int
handle_input_from_worker_undecided(SESSION *sess, WORKER *worker)
{
    if (queue_num_filled(worker->inbuf) >= 1)
        switch (worker->inbuf->data[worker->inbuf->front]) {
        case '@': /* FASTQ */
            worker->input_handler = handle_input_from_worker_fastq;
            return handle_input_from_worker_fastq(sess, worker);
        case '>': /* FASTA */
            worker->input_handler = handle_input_from_worker_fasta;
            return handle_input_from_worker_fasta(sess, worker);
        default:
            error(1, "Unknown output format from worker: the first letter "
                     "is not '@' or '>'.\n");
        }

    return 0;
}

static int
launch_workers(SESSION *sess)
{
    int i, j;

    sess->workers = malloc(sizeof(WORKER) * sess->num_workers);
    if (sess->workers == NULL)
        return -1;

    memset(sess->workers, 0, sizeof(WORKER) * sess->num_workers);

    for (i = 0; i < sess->num_workers; i++) {
        int stdin_pipes[2], stdout_pipes[2];

        sess->workers[i].workerid = i;

        if (pipe(stdin_pipes) != 0 || pipe(stdout_pipes) != 0)
            error(1, "Failed to create new pipes.\n");

        if ((sess->workers[i].pid = fork()) == 0) {
            // child
            free(sess->workers);

            dup2(stdin_pipes[0], STDIN_FILENO);
            dup2(stdout_pipes[1], STDOUT_FILENO);
            close(stdin_pipes[1]);
            close(stdout_pipes[0]);

            for (j = 0; j < i; j++)
                if (sess->workers[j].status == STATUS_RUNNING) {
                    close(sess->workers[j].stdin_fd);
                    close(sess->workers[j].stdout_fd);
                }

            if (sess->command != NULL)
                execl("/bin/sh", "sh", "-c", sess->command, (char *)NULL);
            else
                execvp(sess->args[0], sess->args);
            error(1, "Failed to invoke a worker process.\n");
        }

        sess->workers[i].stdin_fd = stdin_pipes[1];
        sess->workers[i].stdout_fd = stdout_pipes[0];
        close(stdin_pipes[0]);
        close(stdout_pipes[1]);

        sess->workers[i].input_handler = handle_input_from_worker_undecided;
        sess->workers[i].status = STATUS_RUNNING;
        sess->running_workers++;
    }

    for (i = 0; i < sess->num_workers; i++)
        if (sess->workers[i].status == STATUS_RUNNING) {
            set_io_nonblocking(sess->workers[i].stdin_fd);
            set_io_nonblocking(sess->workers[i].stdout_fd);
        }

    return 0;
}

static int
main_loop(SESSION *sess)
{
    fd_set rfds, wfds, exfds;
    int i, stdin_closed;

    sess->inbuf = queue_new(INLET_QUEUE_SIZE);
    sess->outbuf = queue_new(OUTLET_QUEUE_SIZE);
    if (sess->inbuf == NULL || sess->outbuf == NULL)
        error(1, "Can't allocate the mainstream queues.\n");

    for (i = 0; i < sess->num_workers; i++) {
        sess->workers[i].inbuf = queue_new(INQUEUE_SIZE);
        sess->workers[i].outbuf = queue_new(OUTQUEUE_SIZE);
        if (sess->workers[i].inbuf == NULL || sess->workers[i].outbuf == NULL)
            error(1, "Failed to allocate buffer queues.\n");
    }

    FD_ZERO(&rfds);
    FD_ZERO(&wfds);
    FD_ZERO(&exfds);

    stdin_closed = 0;

    for (;;) {
        int maxfd=-1;
        struct timeval tv, *timeout;

        if (stdin_closed || is_queue_full(sess->inbuf))
            FD_CLR(STDIN_FILENO, &rfds);
        else {
            FD_SET(STDIN_FILENO, &rfds);
            maxfd = FD_MAX(maxfd, STDIN_FILENO);
        }

        if (!is_queue_empty(sess->outbuf)) {
            FD_SET(STDOUT_FILENO, &wfds);
            maxfd = FD_MAX(maxfd, STDOUT_FILENO);
        }
        else
            FD_CLR(STDOUT_FILENO, &wfds);

        for (i = 0; i < sess->num_workers; i++) {
            WORKER *w=&sess->workers[i];

            if (!IS_STDOUT_ALIVE(w) || is_queue_full(w->inbuf))
                FD_CLR(w->stdout_fd, &rfds);
            else {
                FD_SET(w->stdout_fd, &rfds);
                maxfd = FD_MAX(maxfd, w->stdout_fd);
            }

            if (IS_STDIN_ALIVE(w) && !is_queue_empty(w->outbuf)) {
                FD_SET(w->stdin_fd, &wfds);
                maxfd = FD_MAX(maxfd, w->stdin_fd);
            }
            else
                FD_CLR(w->stdin_fd, &wfds);
        }

        if (maxfd < 0) {
            if (sess->running_workers <= 0)
                break;

            tv.tv_sec = 0;
            tv.tv_usec = 50000;
            timeout = &tv;
        }
        else
            timeout = NULL;

        if (select(maxfd + 1, &rfds, &wfds, &exfds, timeout) == -1) {
            if (errno == EINTR)
                continue;

            perror("echidna");
            error(1, "Error on select()\n");
        }

        if (FD_ISSET(STDIN_FILENO, &rfds)) {
            char *bufstart;
            qsize_t bufsize;
            ssize_t rsize;

            bufsize = queue_num_continuous_vacant(sess->inbuf, &bufstart);

            if ((rsize = read(STDIN_FILENO, bufstart, bufsize)) < 0)
                error(1, "Error on reading from stdin.\n");

            if (rsize == 0)
                stdin_closed = 1;
            else
                queue_queued(sess->inbuf, rsize);

            sess->input_handler(sess);
        }

        if (FD_ISSET(STDOUT_FILENO, &wfds)) {
            char *bufstart;
            qsize_t bufsize;
            ssize_t wsize;

            bufsize = queue_num_continuous_filled(sess->outbuf, &bufstart);

            if ((wsize = write(STDOUT_FILENO, bufstart, bufsize)) < 0) {
                if (errno != EAGAIN)
                    error(1, "Error on writing to stdout.\n");
            }

            queue_consumed(sess->outbuf, wsize);

            for (i = 0; i < sess->num_workers; i++) {
                WORKER *w=&sess->workers[i];
                if (!is_queue_empty(w->inbuf))
                    w->input_handler(sess, w);
            }
        }

        for (i = 0; i < sess->num_workers; i++) {
            WORKER *w=&sess->workers[i];

            if (FD_ISSET(w->stdout_fd, &rfds)) {
                char *bufstart;
                qsize_t bufsize;
                ssize_t rsize;

                bufsize = queue_num_continuous_vacant(w->inbuf, &bufstart);

                if ((rsize = read(w->stdout_fd, bufstart, bufsize)) < 0)
                    error(1, "Error on reading from worker %d.\n", i);

                if (rsize == 0)
                    w->status = STATUS_TERMINATED;
                else
                    queue_queued(w->inbuf, rsize);

                w->input_handler(sess, w);
            }

            if (FD_ISSET(w->stdin_fd, &wfds)) {
                char *bufstart;
                qsize_t bufsize;
                ssize_t wsize;

                bufsize = queue_num_continuous_filled(w->outbuf, &bufstart);

                if ((wsize = write(w->stdin_fd, bufstart, bufsize)) < 0) {
                    if (errno != EAGAIN)
                        error(1, "Error on writing to worker %d.\n", i);
                }

                queue_consumed(w->outbuf, wsize);

                sess->input_handler(sess);
            }
            else if (stdin_closed && w->status == STATUS_RUNNING &&
                     is_queue_empty(w->outbuf)) {
                w->status = STATUS_FLUSHING;
                if (close(w->stdin_fd) == -1)
                    perror("close");
            }
        }
    }

    for (i = 0; i < sess->num_workers; i++) {
        queue_destroy(sess->workers[i].inbuf);
        queue_destroy(sess->workers[i].outbuf);
    }

    queue_destroy(sess->inbuf);
    queue_destroy(sess->outbuf);

    return 0;
}

static void
usage(char *command)
{
    printf("\
Usage: %s [options] [command]\n\
\n\
Options:\n\
  -p, --processes=n     invoke n worker processes (default 4)\n\
  -c, --command=\"cmd\"   invoke a shell command\n\
  -h, --help            display this help\n\
\n\
Report bugs to Hyeshik Chang <hyeshik@snu.ac.kr>\n", command);
}

int
main(int argc, char **argv)
{
    SESSION session;
    int c;

    memset(&session, 0, sizeof(SESSION));
    session.num_workers = DEFAULT_NUM_WORKERS;
    session.input_handler = handle_input_from_stdin_undecided;
    global_session = &session;

    for (;;) {
        static struct option long_options[]={
            {"command",     required_argument, 0, 'c'},
            {"help",        no_argument,       0, 'h'},
            {"processes",   required_argument, 0, 'p'},
            {0, 0, 0, 0}
        };
        int option_index=0;
  
        c = getopt_long(argc, argv, "c:hp:", long_options, &option_index);
  
        if (c == -1)
            break;
  
        switch (c) {
        case 'p':
            session.num_workers = atoi(optarg);
            if (session.num_workers < 1)
                error(1, "invalid thread number: %s\n", optarg);
            break;
        case 'c':
            session.command = optarg;
            break;
        case 'h':
        default:
            usage(argv[0]);
            return 0;
        }
    }

    if (session.command == NULL) {
        if (optind >= argc) {
            error(0, "command is not supplied.\n\n");
            usage(argv[0]);
            return 0;
        }

        session.args = argv + optind;
    }

    signal(SIGCHLD, sigchld_hdl);

    launch_workers(&session);

    set_io_nonblocking(STDIN_FILENO);
    set_io_nonblocking(STDOUT_FILENO);

    main_loop(&session);

    signal(SIGCHLD, SIG_DFL);

    return 0;
}
