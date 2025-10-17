use std::fmt::Display;

use crate::fold::{AtFaceCorner, AtHalfEdgeField, Edge, Face, FaceCorner, Frame, FrameAttribute, HalfEdge, Vertex};

/// An error that happens when checking if a frame is a manifold
#[derive(Clone, Debug)]
#[cfg_attr(test, derive(PartialEq))] // Really no point in this derivation outside of tests
pub enum ManifoldError {
    /// An edge is adjacent to at least 3 faces, which is too many.
    EdgeAdjacentToTooManyFaceCorners { edge: Edge, face_corners: [Vec<FaceCorner>; 2] },
    /// A vertex's face corners cannot be sequenced in a fan/loop
    VertexFacesDoNotFormFanOrLoop { vertex: Vertex, face_corners: Vec<FaceCorner> }
}

impl Display for ManifoldError {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        write!(f, "manifold error: ")?;
        match self {
            Self::EdgeAdjacentToTooManyFaceCorners { edge, face_corners: faces } => {
                write!(f, "edge {edge} is adjacent to 3 or more faces ({faces:?})")
            }
            Self::VertexFacesDoNotFormFanOrLoop { vertex, face_corners: faces } => {
                write!(f, "vertex {vertex}'s faces ({faces:?}) do not form a fan or loop")
            }
        }
    }
}

/// An error that happens when checking if a frame is an orientable manifold
//#[derive(Clone, Debug)]
//#[cfg_attr(test, derive(PartialEq))] // Really no point in this derivation outside of tests
//pub enum OrientableError {
//    /// An half-edge is adjacent to at least 2 faces, which is too many.
//    HalfEdgeAdjacentToTooManyFaces { half_edge: HalfEdge, faces: Vec<Face> },
//    /// A vertex's faces cannot be sequenced in a fan/loop
//    VertexFacesDoNotFormFanOrLoop { vertex: Vertex, faces: Vec<Face> }
//}
//
//impl Display for OrientableError {
//    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
//        write!(f, "manifold error: ")?;
//        match self {
//            Self::HalfEdgeAdjacentToTooManyFaces { half_edge, faces } => {
//                write!(f, "half-edge {half_edge} is adjacent to 2 or more faces ({faces:?})")
//            }
//            Self::VertexFacesDoNotFormFanOrLoop { vertex, faces } => {
//                write!(f, "vertex {vertex}'s faces ({faces:?}) do not form a fan or loop")
//            }
//        }
//    }
//}

impl Frame {
    /// Tries to convert this frame into a manifold (possibly with boundary) representation.
    /// On success, every vertex will have its half-edges listed cyclically according to the fan/loop of its faces,
    /// and `FrameAttribute::Manifold` is added to the frame attributes. Nothing else changes.
    /// 
    /// In the case that a vertex is adjacent to some edges with 0 faces attached to them,
    /// the corresponding half-edges go to the end of the list in their previous order.
    /// 
    /// Requirements:
    /// * Each edge must be adjacent to at most 2 faces
    /// * In each vertex, the faces must be sequenceable in such a way that adjacent faces share an edge.
    ///     (in other words, the faces form a fan or loop. Double-cones are not allowed.)
    /// 
    /// Note that a frame with no faces is always a manifold according to this definition.
    pub fn try_into_manifold(mut self) -> Result<Self, Vec<ManifoldError>> {
        self.add_attribute_unchecked(FrameAttribute::Manifold);

        // Welp; now to check it.
        let (vh, ev, ef, fh) = if let (Some(vh), Some(ev), Some(ef), Some(fh)) =
            (self.vertices_half_edges.as_mut(), self.edges_vertices.as_ref(), self.edges_face_corners.as_ref(), self.faces_half_edges.as_ref())
        {
            (vh, ev, ef, fh)
        } else { return Ok(self) };

        let mut errors = vec![];
        for (e, faces) in ef.iter_enumerated() {
            if faces[0].len() + faces[1].len() > 2 {
                errors.push(ManifoldError::EdgeAdjacentToTooManyFaceCorners { edge: e, face_corners: faces.clone() })
            }
        }

        for (v, half_edges) in vh.iter_mut_enumerated() {
            // Find a start: an edge adjacent to exactly 1 face.
            // Otherwise, just take an edge adjacent to 2 faces; otherwise skip
            let start = half_edges.iter().position(|&h|
                ef.at(h).len() + ef.at(h.flipped()).len() == 1);
            let start = start.or_else(|| half_edges.iter().position(|&h|
                ef.at(h).len() + ef.at(h.flipped()).len() == 2));
            let start_half_edge = if let Some(s) = start { half_edges[s] } else { continue };

            let mut curr_half_edge = start_half_edge;
            let mut curr_face_corner = None;
            let mut order = vec![];
            let mut unaccounted = half_edges.iter().copied().map(Some).collect::<Vec<_>>();
            // Build the fan
            //let mut iters = 0;
            println!("vertex: {}", v);
            loop {
                //if iters >= 10 {
                //    panic!("Too many iterations!");
                //}
                //iters += 1;
                println!("unaccounted: {:?}, h: {}, c: {:?}", unaccounted, curr_half_edge, curr_face_corner);
                order.push(curr_half_edge);
                unaccounted.iter_mut().find(|h| **h == Some(curr_half_edge)).map(|h| *h = None);
                println!("new unaccounted: {:?}", unaccounted);

                // Next half-edge
                let corner = ef.pair_at(curr_half_edge).iter().enumerate()
                    .flat_map(|(flip, cs)|
                        cs.iter().map(move |&c| (flip != 0, if flip != 0 { c.next(fh) } else { c })))
                    .find(|&(_, c)| curr_face_corner.map(|c2| c != c2).unwrap_or(true));
                let (flip, corner) = if let Some(c) = corner { c } else { break };
                curr_half_edge = if flip {
                    fh.at(corner)
                } else {
                    fh.at(corner.prev(fh)).flipped()
                };

                curr_face_corner = Some(corner);
                println!("new h: {}, c: {:?}", curr_half_edge, curr_face_corner);

                if curr_half_edge == start_half_edge { break; }
            }

            // Done. All other edges had better be adjacent to 0 faces.
            let unaccounted = unaccounted.iter().flatten().collect::<Vec<_>>();
            if unaccounted.iter().any(|h| ef[h.edge()].iter().flatten().count() > 0) {
                errors.push(ManifoldError::VertexFacesDoNotFormFanOrLoop {
                    vertex: v,
                    face_corners: half_edges.iter().flat_map(|&h| ef.at(h).iter()).copied().collect::<Vec<_>>()
                });
            }

            *half_edges = order;
            half_edges.extend(unaccounted);
        }

        if errors.is_empty() { Ok(self) } else { Err(errors) }
    }
}

#[cfg(test)]
mod test {
    use typed_index_collections::{ti_vec, TiSlice, TiVec};

    use crate::fold::{Frame, FrameAttribute};
    use crate::fold::{Vertex as V, Edge as E, Face as F, HalfEdge as H, FaceCorner as C};
    use crate::manifold::ManifoldError;

    fn assert_manifold_vh_maybe_isomorphic(mut a: TiVec<V, Vec<H>>, mut b: TiVec<V, Vec<H>>, fan_lens: TiVec<V, (usize, bool)>, isomorphic: bool) {
        // Canonicalize first
        for (hs, &(fan_len, is_loop)) in a.iter_mut().zip(fan_lens.iter()).chain(b.iter_mut().zip(fan_lens.iter())) {
            if is_loop {
                *hs = (0..hs.len()).flat_map(|i| [(i, false), (i, true)])
                    .map(|(i, flip)| {
                        let mut result = hs.clone();
                        result[0..fan_len].rotate_left(i);
                        if flip { result[0..fan_len].reverse(); }
                        result
                    })
                    .min().unwrap();
            } else {
                let mut reversed = hs.clone();
                reversed[0..fan_len].reverse();
                *hs = hs.clone().min(reversed);
            }
        }

        if isomorphic { assert_eq!(a, b) } else { assert_ne!(a, b) };
    }

    fn assert_manifold_vh_isomorphic(a: TiVec<V, Vec<H>>, b: TiVec<V, Vec<H>>, fan_lens: TiVec<V, (usize, bool)>) {
        assert_manifold_vh_maybe_isomorphic(a, b, fan_lens, true);
    }

    fn assert_manifold_vh_not_isomorphic(a: TiVec<V, Vec<H>>, b: TiVec<V, Vec<H>>, fan_lens: TiVec<V, (usize, bool)>) {
        assert_manifold_vh_maybe_isomorphic(a, b, fan_lens, false);
    }

    #[test]
    fn test_try_into_manifold() {
        use ManifoldError::*;

        // Empty frame
        let manifold = Frame {
            ..Default::default()
        }.try_into_manifold();
        let manifold = manifold.unwrap();
        assert_eq!(manifold.frame_attributes, vec![FrameAttribute::Manifold]);
        assert_eq!(manifold.vertices_half_edges, None);

        // No faces. Just a 3-star.
        let manifold = Frame {
            ..Default::default()
        }.with_topology_vh_fh(Some(ti_vec![
            vec![H(0), H(2), H(4)],
            vec![H(1)],
            vec![H(3)],
            vec![H(5)],
        ]), None).try_into_manifold();
        assert_manifold_vh_isomorphic(manifold.unwrap().vertices_half_edges.unwrap(), ti_vec![
            vec![H(0), H(2), H(4)],
            vec![H(1)],
            vec![H(3)],
            vec![H(5)],
        ], ti_vec![
            (0, false),
            (0, false),
            (0, false),
            (0, false),
        ]);

        // A simple triangle.
        let manifold = Frame {
            ..Default::default()
        }.with_topology_vh_fh(Some(ti_vec![
            vec![H(5), H(0)],
            vec![H(1), H(2)],
            vec![H(3), H(4)],
        ]), Some(ti_vec![
            vec![H(0), H(2), H(4)],
        ])).try_into_manifold();
        assert_manifold_vh_isomorphic(manifold.unwrap().vertices_half_edges.unwrap(), ti_vec![
            vec![H(5), H(0)],
            vec![H(1), H(2)],
            vec![H(3), H(4)],
        ], ti_vec![
            (2, false),
            (2, false),
            (2, false),
        ]);

        // A simple triangle, with some edges sticking out.
        let manifold = Frame {
            ..Default::default()
        }.with_topology_vh_fh(Some(ti_vec![
            vec![H(5), H(0)],
            vec![H(1), H(2)],
            vec![H(3), H(4), H(6), H(10), H(8)],
            vec![H(7)],
            vec![H(9)],
            vec![H(11)],
        ]), Some(ti_vec![
            vec![H(0), H(2), H(4)],
        ])).try_into_manifold();
        assert_manifold_vh_isomorphic(manifold.unwrap().vertices_half_edges.unwrap(), ti_vec![
            vec![H(5), H(0)],
            vec![H(1), H(2)],
            vec![H(3), H(4), H(6), H(10), H(8)],
            vec![H(7)],
            vec![H(9)],
            vec![H(11)],
        ], ti_vec![
            (2, false),
            (2, false),
            (2, false),
            (0, false),
            (0, false),
            (0, false),
        ]);

        // A doubled square
        let manifold = Frame {
            ..Default::default()
        }.with_topology_vh_fh(Some(ti_vec![
            vec![H(7), H(0)],
            vec![H(1), H(2)],
            vec![H(3), H(4)],
            vec![H(5), H(6)]
        ]), Some(ti_vec![
            vec![H(0), H(2), H(4), H(6)],
            vec![H(7), H(5), H(3), H(1)],
        ])).try_into_manifold();
        assert_manifold_vh_isomorphic(manifold.unwrap().vertices_half_edges.unwrap(), ti_vec![
            vec![H(7), H(0)],
            vec![H(1), H(2)],
            vec![H(3), H(4)],
            vec![H(5), H(6)],
        ], ti_vec![
            (2, true),
            (2, true),
            (2, true),
            (2, true),
        ]);

        // A square cross. Note that the order of vertices in the apex matters a bit.
        let manifold = Frame {
            ..Default::default()
        }.with_topology_vh_fh(Some(ti_vec![
            vec![H(0), H(7), H(8)],
            vec![H(1), H(2), H(10)],
            vec![H(3), H(4), H(12)],
            vec![H(5), H(6), H(14)],
            vec![H(9), H(13), H(11), H(15)],
        ]), Some(ti_vec![
            vec![H(0), H(10), H(9)],
            vec![H(2), H(12), H(11)],
            vec![H(4), H(14), H(13)],
            vec![H(6), H(8), H(15)],
        ])).try_into_manifold();
        assert_manifold_vh_isomorphic(manifold.clone().unwrap().vertices_half_edges.unwrap(), ti_vec![
            vec![H(7), H(8) , H(0)],
            vec![H(1), H(10), H(2)],
            vec![H(3), H(12), H(4)],
            vec![H(5), H(14), H(6)],
            vec![H(9), H(11), H(13), H(15)],
        ], ti_vec![
            (3, false),
            (3, false),
            (3, false),
            (3, false),
            (4, true),
        ]);
        assert_manifold_vh_not_isomorphic(manifold.unwrap().vertices_half_edges.unwrap(), ti_vec![
            vec![H(7), H(8) , H(0)],
            vec![H(1), H(10), H(2)],
            vec![H(3), H(12), H(4)],
            vec![H(5), H(14), H(6)],
            vec![H(9), H(13), H(11), H(15)], // not allowed
        ], ti_vec![
            (3, false),
            (3, false),
            (3, false),
            (3, false),
            (4, true),
        ]);

        // A slitted face
        // 3----2----2
        // |         |
        // |  7-6-6  |
        // 3  7   5  1
        // |  4-4-5  |
        // |,8       |
        // 0----0----1
        let manifold = Frame {
            ..Default::default()
        }.with_topology_vh_fh(Some(ti_vec![
            vec![H(0), H(7), H(16)],
            vec![H(1), H(2)],
            vec![H(3), H(4)],
            vec![H(5), H(6)],
            vec![H(8), H(15), H(17)],
            vec![H(9), H(10)],
            vec![H(11), H(12)],
            vec![H(13), H(14)],
        ]), Some(ti_vec![
            vec![H(0), H(2), H(4), H(6), H(16), H(15), H(13), H(11), H(9), H(17)],
            vec![H(8), H(10), H(12), H(14)],
        ])).try_into_manifold();
        assert_manifold_vh_isomorphic(manifold.unwrap().vertices_half_edges.unwrap(), ti_vec![
            vec![H(0), H(16), H(7)],
            vec![H(1), H(2)],
            vec![H(3), H(4)],
            vec![H(5), H(6)],
            vec![H(8), H(15), H(17)],
            vec![H(9), H(10)],
            vec![H(11), H(12)],
            vec![H(13), H(14)],
        ], ti_vec![
            (3, false),
            (2, false),
            (2, false),
            (2, false),
            (3, true),
            (2, true),
            (2, true),
            (2, true),
        ]);

        // A möbius face. Noticable, a face duplicates a half-edge.
        // Edges go 0->1->2->3->0 and 0->2.
        // Believe it or not, it's a manifold.
        // »---0       3--3»
        //     |"0   ."     
        //     4  "."       
        //     | 2" ".      
        // >---2"     "1--1>  
        let manifold = Frame {
            ..Default::default()
        }.with_topology_vh_fh(Some(ti_vec![
            vec![H(7), H(0), H(8)],
            vec![H(1), H(2)],
            vec![H(3), H(4), H(9)],
            vec![H(5), H(6)],
        ]), Some(ti_vec![
            vec![H(0), H(2), H(9), H(7), H(5), H(9)],
            vec![H(0), H(2), H(4), H(6)],
        ])).try_into_manifold();
        assert_manifold_vh_isomorphic(manifold.unwrap().vertices_half_edges.unwrap(), ti_vec![
            vec![H(7), H(0), H(8)],
            vec![H(1), H(2)],
            vec![H(3), H(4), H(9)],
            vec![H(5), H(6)],
        ], ti_vec![
            (3, true),
            (2, true),
            (3, true),
            (2, true),
        ]);

        // A möbius strip. This time, no face by itself violates orientability, but the whole does. Also, there's a boundary.
        // Edges go 0->1->2->3->0 and 0->2 and 3->1.
        // Believe it or not, it's a manifold.
        // »---0       3--3»
        //     |"0   ."|    
        //     4  "."  5    
        //     | 2" ". |    
        // >---2"     "1--1>  
        let manifold = Frame {
            ..Default::default()
        }.with_topology_vh_fh(Some(ti_vec![
            vec![H(7), H(0), H(8)],
            vec![H(1), H(2), H(11)],
            vec![H(3), H(4), H(9)],
            vec![H(5), H(6), H(10)],
        ]), Some(ti_vec![
            vec![H(0), H(11), H(5), H(9)],
            vec![H(6), H(8), H(3), H(11)],
        ])).try_into_manifold();
        assert_manifold_vh_isomorphic(manifold.unwrap().vertices_half_edges.unwrap(), ti_vec![
            vec![H(7), H(8) , H(0)],
            vec![H(1), H(11), H(2)],
            vec![H(3), H(9) , H(4)],
            vec![H(5), H(10), H(6)],
        ], ti_vec![
            (3, false),
            (3, false),
            (3, false),
            (3, false),
        ]);

        // A doubled-up edge
        let manifold = Frame {
            ..Default::default()
        }.with_topology_vh_fh(Some(ti_vec![
            vec![H(0), H(7), H(8), H(10)],
            vec![H(1), H(2)],
            vec![H(3), H(4), H(9), H(11)],
            vec![H(5), H(6)],
        ]), Some(ti_vec![
            vec![H(0), H(2), H(9)],
            vec![H(4), H(6), H(10)],
            vec![H(8), H(11)],
        ])).try_into_manifold();
        assert_manifold_vh_isomorphic(manifold.unwrap().vertices_half_edges.unwrap(), ti_vec![
            vec![H(0), H(8), H(10), H(7)],
            vec![H(1), H(2)],
            vec![H(3), H(9), H(11), H(4)],
            vec![H(5), H(6)],
        ], ti_vec![
            (4, false),
            (2, false),
            (4, false),
            (2, false),
        ]);

        // Three triangles from a single edge (not a manifold)
        let not_a_manifold = Frame {
            ..Default::default()
        }.with_topology_vh_fh(Some(ti_vec![
            vec![H(0), H(5), H(9), H(13)],
            vec![H(1), H(2), H(6), H(10)],
            vec![H(3), H(4)],
            vec![H(7), H(8)],
            vec![H(11), H(12)],
        ]), Some(ti_vec![
            vec![H(0), H(2), H(4)],
            vec![H(0), H(6), H(8)],
            vec![H(0), H(10), H(12)],
        ])).try_into_manifold();
        assert!(matches!(not_a_manifold.unwrap_err()[0], EdgeAdjacentToTooManyFaceCorners { edge: E(0), face_corners: _ }));

        // A square cross, but two opposite faces are missing. Not a manifold.
        let not_a_manifold = Frame {
            ..Default::default()
        }.with_topology_vh_fh(Some(ti_vec![
            vec![H(0), H(7), H(8)],
            vec![H(1), H(2), H(10)],
            vec![H(3), H(4), H(12)],
            vec![H(5), H(6), H(14)],
            vec![H(9), H(13), H(11), H(15)],
        ]), Some(ti_vec![
            vec![H(2), H(12), H(11)],
            vec![H(6), H(8), H(15)],
        ])).try_into_manifold();
        assert!(matches!(not_a_manifold.unwrap_err()[0], VertexFacesDoNotFormFanOrLoop { vertex: V(4), face_corners: _ }));
    }
}