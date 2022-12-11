use std::collections::HashMap;
use std::*;

use bairstow::rootfinding::{delta, delta2, suppress, suppress_old};
use bairstow::vector2::Vector2;
use pytest::approx;
fn test_suppress1() {
    let vri = Vector2(-2, 0);
    let vrj = Vector2(4, -5);
    let mut vA = Vector2(3, 3);
    let mut vA1 = Vector2(1, 2);
    suppress_old(vA, vA1, vri, vrj);
    let dr_old = delta(vA, vri, Vector2(vA1._x, -(vA1._y)));
    vA = Vector2(3, 3);
    vA1 = Vector2(1, 2);
    suppress(vA, vA1, vri, vrj);
    let dr_new = delta2(vA, vri, vA1);
    assert!(dr_new._x == approx(dr_old._x));
    assert!(dr_new._y == approx(dr_old._y));
}
fn test_suppress2() {
    let vr0 = Vector2(-2, 0);
    let vr1 = Vector2(4, -5);
    let vr2 = Vector2(-1, 3);
    let mut vA = Vector2(3, 3);
    let mut vA1 = Vector2(1, 2);
    suppress_old(vA, vA1, vr0, vr1);
    suppress_old(vA, vA1, vr0, vr2);
    let dr_old = delta(vA, vr0, Vector2(vA1._x, -(vA1._y)));
    vA = Vector2(3, 3);
    vA1 = Vector2(1, 2);
    suppress_old(vA, vA1, vr0, vr2);
    suppress_old(vA, vA1, vr0, vr1);
    let dr_old2 = delta(vA, vr0, Vector2(vA1._x, -(vA1._y)));
    assert!(dr_old2._x == approx(dr_old._x));
    assert!(dr_old2._y == approx(dr_old._y));
    vA = Vector2(3, 3);
    vA1 = Vector2(1, 2);
    suppress(vA, vA1, vr0, vr1);
    suppress(vA, vA1, vr0, vr2);
    let dr_new = delta2(vA, vr0, vA1);
    vA = Vector2(3, 3);
    vA1 = Vector2(1, 2);
    suppress(vA, vA1, vr0, vr2);
    suppress(vA, vA1, vr0, vr1);
    let dr_new2 = delta2(vA, vr0, vA1);
    assert!(dr_new._x == approx(dr_new2._x));
    assert!(dr_new._y == approx(dr_new2._y));
    assert!(dr_new._x == approx(dr_old._x));
    assert!(dr_new._y == approx(dr_old._y));
}
