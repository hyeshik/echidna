/**
 * bufqueue.h
 *  - Circular queue byte string buffer implementation for echidna
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

typedef uint32_t qsize_t;

typedef struct {
    qsize_t front;
    qsize_t rear;
    qsize_t size;
    char data[1];
} QUEUE;


static QUEUE *
queue_new(qsize_t size)
{
    QUEUE *q;

    q = malloc(sizeof(qsize_t)*3 + size);
    if (q == NULL)
        return NULL;

    q->front = q->rear = 0;
    q->size = size;

    return q;
}

static void
queue_destroy(QUEUE *q)
{
    free(q);
}

static inline qsize_t
queue_num_filled(QUEUE *q)
{
    return (q->rear + q->size - q->front) % q->size;
}

static inline qsize_t
queue_num_vacant(QUEUE *q)
{
    return (q->front + q->size - q->rear - 1) % q->size;
}

static inline qsize_t
queue_num_continuous_vacant(QUEUE *q, char **bufstart)
{
    *bufstart = q->data + q->rear;

    if (q->rear >= q->front)
        if (q->front == 0)
            return q->size - q->rear - 1;
        else
            return q->size - q->rear;
    else
        return q->front - q->rear - 1;
}

static inline qsize_t
queue_num_continuous_filled(QUEUE *q, char **bufstart)
{
    *bufstart = q->data + q->front;

    if (q->rear >= q->front)
        return q->rear - q->front;
    else
        return q->size - q->front;
}

static inline int
is_queue_empty(QUEUE *q)
{
    return q->front == q->rear;
}

static inline int
is_queue_full(QUEUE *q)
{
    return (q->rear + 1 + q->size) % q->size == q->front;
}

static inline void
queue_consumed(QUEUE *q, qsize_t step)
{
    q->front = (q->front + step) % q->size;
}

static inline void
queue_queued(QUEUE *q, qsize_t step)
{
    q->rear = (q->rear + step) % q->size;
}

static inline int
queue_put(QUEUE *q, char *data, qsize_t datalen)
{
    qsize_t right_vacant;

    if (datalen > queue_num_vacant(q))
    	return -1;

    right_vacant = q->size - q->rear;
    if (right_vacant >= datalen) /* enough space on the right side */
        memcpy(q->data + q->rear, data, datalen);
    else {
        memcpy(q->data + q->rear, data, right_vacant);
        memcpy(q->data, data + right_vacant, datalen - right_vacant);
    }

    queue_queued(q, datalen);

    return 0;
}

static inline int
queue_transfer(QUEUE *dst, QUEUE *src, qsize_t len)
{
    qsize_t right_contig;
    int ret;

    right_contig = src->size - src->front;
    if (right_contig >= len) /* not fragmented aside */
        ret = queue_put(dst, src->data + src->front, len);
    else if (len <= queue_num_vacant(dst)) {
        if (queue_put(dst, src->data + src->front, right_contig) != 0 ||
                queue_put(dst, src->data, len - right_contig) != 0)
            errx(0, "A queue operation is out of integrity.\n");
        ret = 0;
    }
    else
        ret = -1;

    if (ret == 0)
        queue_consumed(src, len);

    return ret;
}
