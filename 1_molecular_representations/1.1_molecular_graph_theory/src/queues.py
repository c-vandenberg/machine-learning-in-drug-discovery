from queue import Queue
from typing import Any


class FifoQueue(Queue):
    """
    A simple FIFO (first-in, first-out) queue. No changes made to Queue base class, simply created to get author
    familiar with Python Queue classes

    Methods
    -------
    enqueue(element: Any)
        Append element to back of queue.
    dequeue() -> Any
        Remove and return element at front of queue.
    is_empty() -> bool
        Return True if the queue is empty, False otherwise.
    is_full() -> bool
        Return True if the queue is full, False otherwise.
    queue_size() -> int
        Return the approximate size of the queue.
    """

    def enqueue(self, element: Any):
        """
        Append element to back of queue.

        Parameters
        ----------
        element : Any
            The element to be added to the queue.
        """
        self.put(element)

    def dequeue(self) -> Any:
        """
        Remove and return element at front of queue.

        Returns
        -------
        Any
            The element at the front of the queue.
        """
        return self.get()

    def is_empty(self) -> bool:
        """
        Return True if the queue is empty, False otherwise.

        Returns
        -------
        bool
            True if the queue is empty, False otherwise.
        """
        return self.empty()

    def is_full(self) -> bool:
        """
        Return True if the queue is full, False otherwise.

        Returns
        -------
        bool
            True if the queue is full, False otherwise.
        """
        return self.full()

    def queue_size(self) -> int:
        """
        Return the approximate size of the queue.

        Returns
        -------
        int
            The approximate size of the queue.
        """
        return self.qsize()
