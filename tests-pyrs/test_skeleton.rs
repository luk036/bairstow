use std::collections::HashMap;
use std::*;

use bairstow::skeleton::{fib, main};
const __author__: _ = "Wai-Shing Luk";
const __copyright__: _ = "Wai-Shing Luk";
const __license__: _ = "MIT";
fn test_fib() {
    "API Tests";
    assert!(fib(1) == 1);
    assert!(fib(2) == 1);
    assert!(fib(7) == 13);
    // with!(pytest.raises(AssertionError)) //unsupported
    {
        fib(-10);
    }
}
fn test_main<T0>(capsys: T0) {
    "CLI Tests";
    main(vec!["7"]);
    let captured = capsys.readouterr();
    assert!(captured
        .out
        .iter()
        .any(|&x| x == "The 7-th Fibonacci number is 13"));
}
