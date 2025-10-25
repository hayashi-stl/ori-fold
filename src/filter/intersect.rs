use nalgebra::DMatrix;

use crate::{filter::Coordinate, geom::NumEx, Frame};

// Pseudocode for line sweep (documenting for reference)
//
// Transform all coordinates (x, y) into (x + εy, y) for infinitesimal ε to remove vertical lines and incidental same-time events
// Delete all length-0 segments (those where x coordinates are the same now)
// PQ <- new priority queue
// SL <- new AVL tree
// add start and end events of all segments into PQ
// while e <- pop PQ
//     E <- pop all events happening at the same time as e from PQ
//     add e to E
//     SL_E <- remove all segments involved in E from SL
//     for all End(s) in E
//         remove s from SL_E
//     reverse order of SL_E
//     split segments in SL_E
//     if SL_E is empty and E has a Start event Start(s):
//         attempt to insert s into SL
//         remove s from SL
//         if there was a tie by position (don't care about segment angle here)
//             t <- offending segment
//             remove t from SL
//             add t to SL_E
//             split t
//     for all Start(s) in E
//         attempt to insert s into SL_E
//         if there's a tie by angle
//             remove s from SL_E
//             t <- offending segment
//             if End(s) < End(t)
//                 add Intersect(t, t) at same time as End(s) to PQ
//                 delete s
//             else if End(s) = End(t)
//                 delete s
//             else
//                 add Start(s) at same time as End(t) to PQ
//     insert SL_E into SL
//     check intersection between start of SL_E and previous in SL, and add Intersection event to PQ if there's one
//     if SL_E is not empty
//         check intersection between end of SL_E and next in SL, and add Intersection event to PQ if there's one

impl Frame {
    /// Adds vertices at all points of intersections between the edges of `self`.
    /// 
    /// This requires coordinates to be 2D, and will remove all face information. Beware. 
    /// 
    /// In particular:
    /// * Duplicate edges are merged. Even ones that showed up during processing.
    /// * If two vertices have the same coordinates, they get merged. Even ones that showed up during processing.
    /// * If two edges are collinear and intersect in their interiors, each edge is split
    ///     wherever it intersects a boundary point of the other edge. Then duplicate edges are merged.
    /// * If an edge intersects a boundary point of another edge, the first edge is split at that boundary point.
    /// * If two non-collinear edges intersect in their interiors, they're both split at the point of intersection.
    pub fn intersect_all_edges_generic<T: NumEx + Coordinate>(&mut self) {
    }
}