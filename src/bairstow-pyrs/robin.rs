use std::collections::HashMap;
use std::*;

struct SlNode {
    next: ST0,
    data: ST1,
}

impl SlNode {
    fn __init__<T0>(&self, data: T0) {
        // "initialization
        //
        //         Keyword Arguments:
        //             data (type):  description
        //         ";
        self.next = self;
        self.data = data;
    }
}
struct Robin {
    cycle: ST0,
}

impl Robin {
    // "Round Robin
    //
    //     Raises:
    //         StopIteration:  description
    //
    //     Returns:
    //         dtype:  description
    //     ";
    const __slots__: _ = "cycle";
    fn __init__(&self, num_parts: i32) {
        self.cycle = (0..num_parts)
            .iter()
            .map(|k| SlNode(k))
            .collect::<Vec<_>>()
            .collect::<Vec<_>>();
        let mut sl2 = self.cycle[-1];
        for sl1 in self.cycle {
            sl2.next = sl1;
            sl2 = sl1;
        }
    }
    fn exclude<RT>(&self, from_part: i32) -> RT {
        // "iterator
        //
        //         Returns:
        //             robin_iterator
        //         ";
        return robin_iterator(self, from_part);
    }
}
struct robin_iterator {
    cur: ST0,
}

impl robin_iterator {
    fn __init__<T0>(&self, Robin: T0, from_part: i32) {
        // "[summary]
        //
        // Arguments:
        //     Robin (type):  description
        // ";
        self.cur = Robin::cycle[from_part];
    }
    fn __iter__<RT>(&self) -> RT {
        // "iterable
        //
        // Returns:
        //     robin_iterator:  itself
        // ";
        return self;
    }
    fn next<RT>(&self) -> RT {
        // "next
        //
        // Raises:
        //     StopIteration:  description
        //
        // Returns:
        //     robinink:  description
        // ";
        self.cur = self.cur.next;
        if self.cur != self.stop {
            return self.cur.data;
        } else {
            raise!(StopIteration); //unsupported
        }
    }
    fn __next__<RT>(&self) -> RT {
        // "[summary]
        //
        // Returns:
        //     dtype:  description
        // ";
        return self.next();
    }
}
fn main() {
    let r = Robin(5);
    for k in r.exclude(3) {
        println!("{:?} ", k);
    }
}
