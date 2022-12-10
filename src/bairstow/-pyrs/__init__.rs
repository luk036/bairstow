use std::*;
use std::collections::HashMap;

if sys.version_info[..2] >= (3, 8) {
use importlib::metadata::{PackageNotFoundError, version};
} else {
use importlib_metadata::{PackageNotFoundError, version};
}
let try_dummy = { //unsupported
const dist_name: _ = __name__;
const __version__: _ = version(dist_name);
};
let except!(PackageNotFoundError) = { //unsupported
__version__ = "unknown";
};