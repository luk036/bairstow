from typing import List


class SlNode:
    next: "SlNode"
    data: int

    def __init__(self, data: int):
        """_summary_

        Args:
            data (int): _description_
        """
        self.next = self
        self.data = data


class RobinIterator:
    __slots__ = ("cur", "stop")
    cur: SlNode
    stop: SlNode

    def __init__(self, node: SlNode) -> None:
        """_summary_

        Args:
            node (SlNode): _description_
        """
        self.cur = self.stop = node

    def __iter__(self) -> "RobinIterator":
        """_summary_

        Returns:
            RobinIterator: _description_
        """
        return self

    def next(self) -> int:
        """_summary_

        Raises:
            StopIteration: _description_

        Returns:
            int: _description_
        """
        self.cur = self.cur.next
        if self.cur != self.stop:
            return self.cur.data
        else:
            raise StopIteration()

    def __next__(self):
        """_summary_

        Returns:
            _type_: _description_
        """
        return self.next()


class Robin:
    __slots__ = "cycle"
    cycle: List[SlNode]

    def __init__(self, num_parts: int):
        """_summary_

        Args:
            num_parts (int): _description_
        """
        self.cycle = list(SlNode(k) for k in range(num_parts))
        sl2 = self.cycle[-1]
        for sl1 in self.cycle:
            sl2.next = sl1
            sl2 = sl1

    def exclude(self, from_part: int) -> RobinIterator:
        """_summary_

        Args:
            from_part (int): _description_

        Returns:
            RobinIterator: _description_
        """
        return RobinIterator(self.cycle[from_part])


if __name__ == "__main__":
    r = Robin(5)
    for k in r.exclude(3):
        print(k)
