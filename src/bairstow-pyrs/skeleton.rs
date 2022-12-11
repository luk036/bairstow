use std::*;
use std::collections::HashMap;

"
This is a skeleton file that can serve as a starting point for a Python
console script. To run this script uncomment the following lines in the
``[options.entry_points]`` section in ``setup.cfg``::

    console_scripts =
         fibonacci = bairstow.skeleton:run

Then run ``pip install .`` (or ``pip install -e .`` for editable mode)
which will install the command ``fibonacci`` inside your current environment.

Besides console scripts, the header (i.e. until ``_logger``...) of this file can
also be used as template for Python modules.

Note:
    This skeleton file can be safely removed if not needed!

References:
    - https://setuptools.readthedocs.io/en/latest/userguide/entry_point.html
    - https://pip.pypa.io/en/stable/reference/pip_install
";
use bairstow::{__version__};
const __author__: _ = "Wai-Shing Luk";
const __copyright__: _ = "Wai-Shing Luk";
const __license__: _ = "MIT";
const _logger: _ = logging.getLogger(__name__);
fn fib<T0, RT>(n: T0) -> RT {
"Fibonacci example function

    Args:
      n (int): integer

    Returns:
      int: n-th Fibonacci number
    ";
assert!(n > 0);
let (a, b) = (1, 1);
for _i in (0..(n - 1)) {
let (a, b) = (b, (a + b));
}
return a;
}
fn parse_args<T0, RT>(args: T0) -> RT {
"Parse command line parameters

    Args:
      args (List[str]): command line parameters as list of strings
          (for example  ``[\"--help\"]``).

    Returns:
      :obj:`argparse.Namespace`: command line parameters namespace
    ";
let parser = argparse.ArgumentParser("Just a Fibonacci demonstration");
parser.add_argument("--version", "version", "bairstow {ver}".format(__version__));
parser.add_argument("n", "n-th Fibonacci number", int, "INT");
parser.add_argument("-v", "--verbose", "loglevel", "set loglevel to INFO", "store_const", logging.INFO);
parser.add_argument("-vv", "--very-verbose", "loglevel", "set loglevel to DEBUG", "store_const", logging.DEBUG);
return parser.parse_args(args);
}
fn setup_logging<T0>(loglevel: T0)  {
"Setup basic logging

    Args:
      loglevel (int): minimum loglevel for emitting messages
    ";
let logformat = "[%(asctime)s] %(levelname)s:%(name)s:%(message)s";
logging.basicConfig(loglevel, sys.stdout, logformat, "%Y-%m-%d %H:%M:%S");
}
fn main<T0>(args: T0)  {
"Wrapper allowing :func:`fib` to be called with string arguments in a CLI fashion

    Instead of returning the value from :func:`fib`, it prints the result to the
    ``stdout`` in a nicely formatted message.

    Args:
      args (List[str]): command line parameters as list of strings
          (for example  ``[\"--verbose\", \"42\"]``).
    ";
args = parse_args(args);
setup_logging(args.loglevel);
_logger.debug("Starting crazy calculations...");
println!("{:?} ","The {}-th Fibonacci number is {}".format(args.n, fib(args.n)));
_logger.info("Script ends here");
}
fn run()  {
"Calls :func:`main` passing the CLI arguments extracted from :obj:`sys.argv`

    This function can be used as entry point to create console scripts with setuptools.
    ";
main(sys.argv[1..]);
}
fn main() {
run();
}