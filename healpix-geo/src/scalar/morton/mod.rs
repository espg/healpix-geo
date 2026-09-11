//! The `morton` scheme: HEALPix cells as packed-`u64` Morton words.
//!
//! Backed by the [`mortie-core`](https://crates.io/crates/mortie-core) codec.
//! Like `zuniq`, a Morton word encodes its own refinement level, so most
//! functions take no `depth` parameter. Unlike `zuniq`, a raw unsigned sort of
//! Morton words is a Z-order traversal with parents sorting immediately before
//! their children, and containment testing is plain prefix truncation (see
//! [`hierarchy::contains`]).
pub mod conversion;
pub mod hierarchy;
