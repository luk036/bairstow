# Python, Rust, C++

Python: no function overloading
        has default function arguments

```python
def solve(x, tol=1e-8):
    ...

class Partitioner:
    def __init__(G, k = 2):
        ...

def bad_idea(x, lst=[]):  # Don't!
    ...
```

Rust: 

- no default function arguments
- limited function overloading
  through trait
  allow different types but must have same number of arguments

```rust
pub fn cut(beta: f64)

pub fn cut(beta: (f64, Option<f64>))
```



