use std::{hash::Hash, mem};

use exact_number::BasedExpr;
use indexmap::{IndexMap, indexmap};
use nalgebra::{DVector, DefaultAllocator, Dim, Dyn, OVector, Scalar, U1, Vector, VectorView, allocator::Allocator};

use crate::{AtFaceCorner, Edge, EdgeData, EdgesFaceCornersEx, EdgesVerticesEx, FaceCorner, Frame, FrameAttribute, HalfEdge, Vertex, filter::Coordinate, geom::{FloatOrd, NumEx}};

pub trait MergeCoordinate: Sized + Scalar {
    type Hash<'a, D: Dim>: Hash where
            Self: 'a,
            DefaultAllocator: Allocator<D>;

    fn hash<'a, D: Dim>(vector: VectorView<'a, Self, D, U1, Dyn>, epsilon: &Self) -> Self::Hash<'a, D> where
            DefaultAllocator: Allocator<D>;

    /// Finds a vertex in `map` closer than `epsilon` to `key`, or inserts `vertex` into the map and returns that.
    fn at_or_insert<'a, D: Dim>(map: &mut IndexMap<Self::Hash<'a, D>, Vertex>, key: VectorView<'a, Self, D, U1, Dyn>, vertex: Vertex, epsilon: &Self) -> Vertex where
            DefaultAllocator: Allocator<D>;
}

impl MergeCoordinate for f64 {
    type Hash<'a, D: Dim> = OVector<FloatOrd<f64>, D> where
            Self: 'a,
            DefaultAllocator: Allocator<D>;

    fn hash<'a, D: Dim>(vector: VectorView<'a, Self, D, U1, Dyn>, epsilon: &Self) -> Self::Hash<'a, D> where
            DefaultAllocator: Allocator<D>
    {
        if *epsilon == 0.0 {
            vector.map(|c| FloatOrd(c))
        } else {
            vector.map(|c| FloatOrd((c / epsilon).round()))
        }
    }

    fn at_or_insert<'a, D: Dim>(map: &mut IndexMap<Self::Hash<'a, D>, Vertex>, key: VectorView<'a, Self, D, U1, Dyn>, vertex: Vertex, epsilon: &Self) -> Vertex where
                DefaultAllocator: Allocator<D>
    {
        if *epsilon == 0.0 {
            return *map.entry(Self::hash(key, epsilon)).or_insert(vertex);
        }

        todo!()
    }
}

impl MergeCoordinate for BasedExpr {
    type Hash<'a, D: Dim> = VectorView<'a, Self, D, U1, Dyn> where
            Self: 'a,
            DefaultAllocator: Allocator<D>;

    /// `epsilon` is not used for exact coordinates.
    fn hash<'a, D: Dim>(vector: VectorView<'a, Self, D, U1, Dyn>, _epsilon: &Self) -> Self::Hash<'a, D> where
            DefaultAllocator: Allocator<D>
    {
        vector
    }

    /// `epsilon` is not used for exact coordinates.
    fn at_or_insert<'a, D: Dim>(map: &mut IndexMap<Self::Hash<'a, D>, Vertex>, key: VectorView<'a, Self, D, U1, Dyn>, vertex: Vertex, epsilon: &Self) -> Vertex where
                DefaultAllocator: Allocator<D>
    {
        *map.entry(key).or_insert(vertex)
    }
}

pub struct PointMerger<'a, T: 'a + MergeCoordinate, D: Dim> where
        DefaultAllocator: Allocator<D>
{
    points: IndexMap<T::Hash<'a, D>, Vertex>,
    merge_map: IndexMap<Vertex, Vertex>,
}

impl<'a, T: 'a + MergeCoordinate, D: Dim> PointMerger<'a, T, D> where
        DefaultAllocator: Allocator<D>
{
    pub fn new(vertices: impl IntoIterator<Item = Vertex>, mut mapping: impl FnMut(Vertex) -> VectorView<'a, T, D, U1, Dyn>, epsilon: T) -> Self {
        let mut this = Self { points: indexmap!{}, merge_map: indexmap!{} };
        for vertex in vertices {
            let to = T::at_or_insert(&mut this.points, mapping(vertex), vertex, &epsilon);
            this.merge_map.insert(vertex, to);
        }
        this
    }

    /// Extracts the merge map (where each vertex should get merged into)
    pub fn into_merge_map(self) -> IndexMap<Vertex, Vertex> {
        self.merge_map
    }
}

impl Frame {
    /// Splits an edge at point `at`.
    /// If there are no vertex coordinates of type `T`, `at` gets ignored.
    /// The new vertex has the same data as vertex `edges_vertices[edge][0]`,
    /// except for coordinates. The new edge has the same data as `edge`
    /// except for the vertices it's connected to.
    /// 
    /// Returns the new vertex `v` created, and the new edge created from the split,
    /// which goes from `v` to `edges_vertices[edge][1]`.
    /// 
    /// This is a purely topological operation and might not preserve
    /// geometric properties such as planarity if you're not careful.
    /// It does preserve manifold-ness and oriented-ness.
    /// 
    /// Requires edge data to exist.
    pub fn split_edge<T: NumEx + Coordinate>(&mut self, edge: Edge, at: DVector<T>) -> (Vertex, Edge) {
        // Add new vertex
        let new_v = self.add_vertex_like(self.edges_vertices.as_ref().unwrap()[edge][0], at);

        // Find where the edge is stored in the end vertex, so it can be replaced later
        let ev = self.edges_vertices.as_ref().unwrap();
        let vh = self.vertices_half_edges.as_ref().unwrap();
        let vh_1 = &vh[ev[edge][1]];
        let vh_1_pos = vh_1.iter().position(|h| h.edge() == edge).expect("inconsistent vertex/edge topology");

        // Add new edge. This is unchecked, so we must fix the topology ourselves.
        let e_data = EdgeData { vertices: [new_v, ev[edge][1]], ..self.edge_data(edge) };
        let new_e = self.add_edge_unchecked(e_data, false);

        // Replace the half-edge
        let ev = self.edges_vertices.as_mut().unwrap();
        let vh = self.vertices_half_edges.as_mut().unwrap();
        vh[ev[edge][1]][vh_1_pos] = HalfEdge::new(new_e, true);
        
        // Connect the split edges together
        ev[edge][1] = new_v;
        vh[new_v] = vec![HalfEdge::new(edge, true), HalfEdge::new(new_e, false)];

        // If there are no faces, we are done.
        let fh = if let Some(fh) = self.faces_half_edges.as_mut() { fh } else { return (new_v, new_e) };
        let ec = self.edges_face_corners.as_mut().unwrap();
        // Unfortunately, if the mesh is oriented and the split edge is adjacent to only one face,
        // then the order of vh[new_v] matters.
        if ec[edge][0].len() == 1 && ec[edge][1].len() == 0 { vh[new_v].reverse() }
        
        // Update faces to use split edge
        let mut corners = ec[edge].iter().enumerate()
            .flat_map(|(i, cs)| cs.iter().map(move |&c| (c, i)))
            .collect::<Vec<_>>();
        // Go in reverse order so indices don't get messed up
        corners.sort();
        corners.reverse();
        let half_edges = new_e.split();
        for (c, i) in corners {
            let insert_index = c.corner() + (1 - i);
            // Fix face corner numbering
            for j in (insert_index..fh[c.face()].len()).rev() {
                let used_h = fh.at(FaceCorner(c.face(), j));
                let used_c = ec.at_mut(used_h).iter_mut().find(|c| c == &&FaceCorner(c.face(), j))
                    .expect("inconsisistent edge/face topology");
                used_c.1 += 1;
            }
            // And add the new half-edge into the face
            fh[c.face()].insert(insert_index, half_edges[i]);
            ec.at_mut(half_edges[i]).push(FaceCorner(c.face(), insert_index));
        }

        (new_v, new_e)
    }

    pub fn split_edges<T: NumEx + Coordinate>(&mut self, splits: impl IntoIterator<Item = (Edge, DVector<T>)>) {
        let vc = T::vertices_coords(self).as_ref();
        let ev = self.edges_vertices.as_ref().unwrap(); // needs to exist to split edges
        let mut splits = splits.into_iter().collect::<Vec<_>>();

        // Need to sort backwards by parameter along edge
        // to avoid splitting an edge where a new edge should already be
        if let Some(vc) = vc {
            splits.sort_by_key(|(e, at)| {
                let coords = ev[*e].map(|v| vc.column(v.0));
                let diff = &coords[0] - &coords[1];
                T::into_sortable((at - &coords[1]).dot(&diff) / diff.norm_squared())
            });
        }

        // Split, split, split!
        for (e, at) in splits {
            self.split_edge(e, at);
        }
    }

    /// Merges `vertex` into `target`, swap-removing `vertex`.
    /// `vertex` and `target` must be different vertices.
    /// 
    /// This effectively replaces all references to `vertex` with `target`,
    /// and then deletes self-loops.
    /// 
    /// This might not preserve manifold-ness. If it preserves manifoldness,
    /// then it preserves orientability. (Note that even if merging vertices causes an
    /// edge to be doubled, those edges are *not* merged here.)
    /// 
    /// If there's an edge connecting the two vertices, edge indices might not
    /// be preserved. If there's a face that only uses edges connecting the two vertices,
    /// face indices might not be preserved. Instead, faces/edges may get swap-removed,
    /// in order from last to first.
    pub fn merge_two_vertices(&mut self, vertex: Vertex, target: Vertex) {
        if let Some(ev) = self.edges_vertices.as_mut() {
            // All edges using `vertex` must now use `target` instead.
            let vh = self.vertices_half_edges.as_mut().unwrap();
            for &h in &vh[vertex] {
                *ev.at_mut(h)[0] = target;
            }
            let vh_vertex = mem::take(&mut vh[vertex]);
            vh[target].extend(vh_vertex);

            if let Some(fh) = self.faces_half_edges.as_mut() {
                // Dissolve self-loops in faces, so that faces no longer use self-loops
                let ec = self.edges_face_corners.as_mut().unwrap();
                let loop_cs = fh.iter_enumerated()
                    .flat_map(|(f, hs)| hs.iter().copied().enumerate().rev()
                        .filter(|&(_, h)| ev.at(h)[0] == ev.at(h)[1])
                        .map(move |(i, _)| FaceCorner(f, i))
                    ).collect::<Vec<_>>();
                for c in loop_cs {
                    // Remove entry from edges_half_faces to keep topology consistent
                    // for later calls to element removal functions
                    let offending_cs = ec.at_mut(fh.at(c));
                    let pos = offending_cs.iter().position(|&c2| c == c2)
                        .expect("inconsistent edge/face topology");
                    offending_cs.remove(pos);

                    // Fix face corner numbering
                    for j in (c.corner() + 1)..fh[c.face()].len() {
                        let used_h = fh.at(FaceCorner(c.face(), j));
                        let used_c = ec.at_mut(used_h).iter_mut().find(|c| c == &&FaceCorner(c.face(), j))
                            .expect("inconsisistent edge/face topology");
                        used_c.1 -= 1;
                    }
                    // Remove self-loop from face
                    fh[c.face()].remove(c.corner());
                }

                // Remove faces that now have 0 edges
                let mut zero_fs = fh.iter_enumerated().filter_map(|(f, hs)| (hs.len() == 0).then_some(f))
                    .collect::<Vec<_>>();
                zero_fs.reverse();
                for f in zero_fs {
                    self.swap_remove_face(f);
                }
            }

            // Remove self-loop edges, now that no faces are using them
            let ev = self.edges_vertices.as_ref().unwrap();
            let mut loop_es = ev.iter_enumerated().filter_map(|(e, vs)| (vs[0] == vs[1]).then_some(e))
                .collect::<Vec<_>>();
            loop_es.reverse();
            for e in loop_es {
                self.swap_remove_edge(e);
            }

            // Fix orientation of merged-into vertex, unless no longer a manifold
            let vh = self.vertices_half_edges.as_mut().unwrap();
            let hs = vh[target].clone();
            self.fix_manifold_attributes_on_half_edges(hs);
        }

        // Whew, after all that topology fixing, we can *finally* remove the vertex
        self.swap_remove_vertex(dbg!(vertex));
    }
}

#[cfg(test)]
mod test {
    use indexmap::{indexmap, indexset};
    use nalgebra::{DMatrix, DVector, U2, vector};
    use typed_index_collections::{TiVec, ti_vec};

    use crate::filter::split_merge::PointMerger;
    use crate::{Frame, FrameAttribute};
    use crate::{Vertex as V, Edge as E, HalfEdge as H, FaceCorner as C, Face as F};

    #[test]
    fn test_split_edge() {
        // No faces. Just a 3-star.
        let mut frame = Frame {
            vertices_coords_f64: Some(DMatrix::from_vec(2, 4, vec![
                0.0, 0.0,
                1.0, 0.0,
                0.0, 1.0,
                -1.0, 0.0,
            ])),
            frame_attributes: indexset![FrameAttribute::Manifold, FrameAttribute::Orientable],
            ..Default::default()
        }.with_topology_vh_fh(Some(ti_vec![
            vec![H(0), H(2), H(4)],
            vec![H(1)],
            vec![H(3)],
            vec![H(5)],
        ]), None);
        let (new_vertex, new_edge) = frame.split_edge::<f64>(E(0), DVector::from_vec(vec![0.5, 0.0]));
        assert_eq!(new_vertex, V(4));
        assert_eq!(new_edge, E(3));
        assert_eq!(frame.vertices_coords_f64, Some(DMatrix::from_vec(2, 5, vec![
            0.0, 0.0,
            1.0, 0.0,
            0.0, 1.0,
            -1.0, 0.0,
            0.5, 0.0,
        ])));
        assert_eq!(frame.edges_vertices.as_ref().unwrap()[E(0)], [V(0), new_vertex]);
        assert_eq!(frame.edges_vertices.as_ref().unwrap()[new_edge], [new_vertex, V(1)]);
        frame.assert_topology_consistent();

        // A simple triangle.
        let mut frame = Frame {
            frame_attributes: indexset![FrameAttribute::Manifold, FrameAttribute::Orientable],
            ..Default::default()
        }.with_topology_vh_fh(Some(ti_vec![
            vec![H(5), H(0)],
            vec![H(1), H(2)],
            vec![H(3), H(4)],
        ]), Some(ti_vec![
            vec![H(5), H(3), H(1)],
        ]));
        let (new_vertex, new_edge) = frame.split_edge::<f64>(E(1), DVector::from_vec(vec![]));
        assert_eq!(new_vertex, V(3));
        assert_eq!(new_edge, E(3));
        assert_eq!(frame.edges_vertices.as_ref().unwrap()[E(1)], [V(1), new_vertex]);
        assert_eq!(frame.edges_vertices.as_ref().unwrap()[new_edge], [new_vertex, V(2)]);
        assert_eq!(frame.faces_half_edges.as_ref().unwrap()[F(0)], vec![H(5), new_edge.split()[1], H(3), H(1)]);
        frame.assert_topology_consistent();

        // Now it's oriented the other way. The half-edge order better be right.
        let mut frame = Frame {
            frame_attributes: indexset![FrameAttribute::Manifold, FrameAttribute::Orientable],
            ..Default::default()
        }.with_topology_vh_fh(Some(ti_vec![
            vec![H(0), H(5)],
            vec![H(2), H(1)],
            vec![H(4), H(3)],
        ]), Some(ti_vec![
            vec![H(0), H(2), H(4)],
        ]));
        let (new_vertex, new_edge) = frame.split_edge::<f64>(E(1), DVector::from_vec(vec![]));
        assert_eq!(new_vertex, V(3));
        assert_eq!(new_edge, E(3));
        assert_eq!(frame.edges_vertices.as_ref().unwrap()[E(1)], [V(1), new_vertex]);
        assert_eq!(frame.edges_vertices.as_ref().unwrap()[new_edge], [new_vertex, V(2)]);
        assert_eq!(frame.faces_half_edges.as_ref().unwrap()[F(0)], vec![H(0), H(2), new_edge.split()[0], H(4)]);
        frame.assert_topology_consistent();

        // A doubled square
        let mut frame = Frame {
            frame_attributes: indexset![FrameAttribute::Manifold, FrameAttribute::Orientable],
            ..Default::default()
        }.with_topology_vh_fh(Some(ti_vec![
            vec![H(7), H(0)],
            vec![H(1), H(2)],
            vec![H(3), H(4)],
            vec![H(5), H(6)]
        ]), Some(ti_vec![
            vec![H(0), H(2), H(4), H(6)],
            vec![H(7), H(5), H(3), H(1)],
        ]));
        let (new_vertex, new_edge) = frame.split_edge::<f64>(E(3), DVector::from_vec(vec![]));
        let [hfwd, hbck] = new_edge.split();
        assert_eq!(new_vertex, V(4));
        assert_eq!(new_edge, E(4));
        assert_eq!(frame.edges_vertices.as_ref().unwrap()[E(3)], [V(3), new_vertex]);
        assert_eq!(frame.edges_vertices.as_ref().unwrap()[new_edge], [new_vertex, V(0)]);
        assert_eq!(frame.faces_half_edges.as_ref().unwrap()[F(0)], vec![H(0), H(2), H(4), H(6), hfwd]);
        assert_eq!(frame.faces_half_edges.as_ref().unwrap()[F(1)], vec![hbck, H(7), H(5), H(3), H(1)]);
        frame.assert_topology_consistent();

        // A slitted face
        // 3----2----2
        // |         |
        // |  7-6-6  |
        // 3  7   5  1
        // |  4-4-5  |
        // |,8       |
        // 0----0----1
        let mut frame = Frame {
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
        ]));
        let (new_vertex, new_edge) = frame.split_edge::<f64>(E(8), DVector::from_vec(vec![]));
        let [hfwd, hbck] = new_edge.split();
        assert_eq!(new_vertex, V(8));
        assert_eq!(new_edge, E(9));
        assert_eq!(frame.edges_vertices.as_ref().unwrap()[E(8)], [V(0), new_vertex]);
        assert_eq!(frame.edges_vertices.as_ref().unwrap()[new_edge], [new_vertex, V(4)]);
        assert_eq!(frame.faces_half_edges.as_ref().unwrap()[F(0)], vec![H(0), H(2), H(4), H(6), H(16), hfwd, H(15), H(13), H(11), H(9), hbck, H(17)]);
        frame.assert_topology_consistent();

        // A möbius face. Noticable, a face duplicates a half-edge.
        // Edges go 0->1->2->3->0 and 0->2.
        // »---0       3--3»
        //     |"0   ."     
        //     4  "."       
        //     | 2" ".      
        // >---2"     "1--1>  
        let mut frame = Frame {
            ..Default::default()
        }.with_topology_vh_fh(Some(ti_vec![
            vec![H(7), H(0), H(8)],
            vec![H(1), H(2)],
            vec![H(3), H(4), H(9)],
            vec![H(5), H(6)],
        ]), Some(ti_vec![
            vec![H(0), H(2), H(9), H(7), H(5), H(9)],
            vec![H(0), H(2), H(4), H(6)],
        ]));
        let (new_vertex, new_edge) = frame.split_edge::<f64>(E(4), DVector::from_vec(vec![]));
        let [hfwd, hbck] = new_edge.split();
        assert_eq!(new_vertex, V(4));
        assert_eq!(new_edge, E(5));
        assert_eq!(frame.edges_vertices.as_ref().unwrap()[E(4)], [V(0), new_vertex]);
        assert_eq!(frame.edges_vertices.as_ref().unwrap()[new_edge], [new_vertex, V(2)]);
        assert_eq!(frame.faces_half_edges.as_ref().unwrap()[F(0)], vec![H(0), H(2), hbck, H(9), H(7), H(5), hbck, H(9)]);
        frame.assert_topology_consistent();

        // Three triangles from a single edge
        let mut frame = Frame {
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
        ]));
        let (new_vertex, new_edge) = frame.split_edge::<f64>(E(0), DVector::from_vec(vec![]));
        let [hfwd, hbck] = new_edge.split();
        assert_eq!(new_vertex, V(5));
        assert_eq!(new_edge, E(7));
        assert_eq!(frame.edges_vertices.as_ref().unwrap()[E(0)], [V(0), new_vertex]);
        assert_eq!(frame.edges_vertices.as_ref().unwrap()[new_edge], [new_vertex, V(1)]);
        assert_eq!(frame.faces_half_edges.as_ref().unwrap()[F(0)], vec![H(0), hfwd, H(2), H(4)]);
        assert_eq!(frame.faces_half_edges.as_ref().unwrap()[F(1)], vec![H(0), hfwd, H(6), H(8)]);
        assert_eq!(frame.faces_half_edges.as_ref().unwrap()[F(2)], vec![H(0), hfwd, H(10), H(12)]);
        frame.assert_topology_consistent();
    }

    #[test]
    fn test_merge_two_vertices() {
        // No faces. Just a 3-star.
        let mut frame = Frame {
            frame_attributes: indexset![FrameAttribute::Manifold, FrameAttribute::Orientable],
            ..Default::default()
        }.with_topology_vh_fh(Some(ti_vec![
            vec![H(0), H(2), H(4)],
            vec![H(1)],
            vec![H(3)],
            vec![H(5)],
        ]), None);
        frame.merge_two_vertices(V(1), V(2));
        assert_eq!(frame.frame_attributes, indexset![FrameAttribute::Manifold, FrameAttribute::Orientable]);
        assert_eq!(frame.edges_vertices, Some(ti_vec![
            [V(0), V(2)],
            [V(0), V(2)],
            [V(0), V(1)],
        ]));
        frame.assert_topology_consistent();

        // A simple square.
        let mut frame = Frame {
            frame_attributes: indexset![FrameAttribute::Manifold, FrameAttribute::Orientable],
            ..Default::default()
        }.with_topology_vh_fh(Some(ti_vec![
            vec![H(7), H(0)],
            vec![H(1), H(2)],
            vec![H(3), H(4)],
            vec![H(5), H(6)],
        ]), Some(ti_vec![
            vec![H(7), H(5), H(3), H(1)],
        ]));
        frame.merge_two_vertices(V(2), V(3));
        assert_eq!(frame.frame_attributes, indexset![FrameAttribute::Manifold, FrameAttribute::Orientable]);
        assert_eq!(frame.edges_vertices, Some(ti_vec![
            [V(0), V(1)],
            [V(1), V(2)],
            [V(2), V(0)],
        ]));
        assert_eq!(frame.faces_half_edges, Some(ti_vec![
            vec![H(5), H(3), H(1)],
        ]));
        frame.assert_topology_consistent();

        // A simple square winded the other way. Make sure to preserve orientable representation
        let mut frame = Frame {
            frame_attributes: indexset![FrameAttribute::Manifold, FrameAttribute::Orientable],
            ..Default::default()
        }.with_topology_vh_fh(Some(ti_vec![
            vec![H(0), H(7)],
            vec![H(2), H(1)],
            vec![H(4), H(3)],
            vec![H(6), H(5)],
        ]), Some(ti_vec![
            vec![H(0), H(2), H(4), H(6)],
        ]));
        frame.merge_two_vertices(V(2), V(3));
        assert_eq!(frame.frame_attributes, indexset![FrameAttribute::Manifold, FrameAttribute::Orientable]);
        assert_eq!(frame.edges_vertices, Some(ti_vec![
            [V(0), V(1)],
            [V(1), V(2)],
            [V(2), V(0)],
        ]));
        assert_eq!(frame.faces_half_edges, Some(ti_vec![
            vec![H(0), H(2), H(4)],
        ]));
        frame.assert_topology_consistent();

        // A slitted face
        // 3----2----2
        // |         |
        // |  7-6-6  |
        // 3  7   5  1
        // |  4-4-5  |
        // |,8       |
        // 0----0----1
        let mut frame = Frame {
            frame_attributes: indexset![FrameAttribute::Manifold, FrameAttribute::Orientable],
            ..Default::default()
        }.with_topology_vh_fh(Some(ti_vec![
            vec![H(0), H(16), H(7)],
            vec![H(2), H(1)],
            vec![H(4), H(3)],
            vec![H(6), H(5)],
            vec![H(8), H(15), H(17)],
            vec![H(9), H(10)],
            vec![H(11), H(12)],
            vec![H(13), H(14)],
        ]), Some(ti_vec![
            vec![H(0), H(2), H(4), H(6), H(16), H(15), H(13), H(11), H(9), H(17)],
            vec![H(8), H(10), H(12), H(14)],
        ]));
        frame.merge_two_vertices(V(4), V(0));
        assert_eq!(frame.frame_attributes, indexset![FrameAttribute::Manifold, FrameAttribute::Orientable]);
        assert_eq!(frame.edges_vertices, Some(ti_vec![
            [V(0), V(1)],
            [V(1), V(2)],
            [V(2), V(3)],
            [V(3), V(0)],
            [V(0), V(5)],
            [V(5), V(6)],
            [V(6), V(4)],
            [V(4), V(0)],
        ]));
        assert_eq!(frame.faces_half_edges, Some(ti_vec![
            vec![H(0), H(2), H(4), H(6), H(15), H(13), H(11), H(9)],
            vec![H(8), H(10), H(12), H(14)],
        ]));
        frame.assert_topology_consistent();

        // A möbius face. Noticable, a face duplicates a half-edge.
        // Edges go 0->1->2->3->0 and 0->2.
        // »---0       3--3»
        //     |"0   ."     
        //     4  "."       
        //     | 2" ".      
        // >---2"     "1--1>  
        let mut frame = Frame {
            ..Default::default()
        }.with_topology_vh_fh(Some(ti_vec![
            vec![H(7), H(0), H(8)],
            vec![H(1), H(2)],
            vec![H(3), H(4), H(9)],
            vec![H(5), H(6)],
        ]), Some(ti_vec![
            vec![H(0), H(2), H(9), H(7), H(5), H(9)],
            vec![H(0), H(2), H(4), H(6)],
        ]));
        frame.merge_two_vertices(V(2), V(0));
        assert_eq!(frame.edges_vertices, Some(ti_vec![
            [V(0), V(1)],
            [V(1), V(0)],
            [V(0), V(2)],
            [V(2), V(0)],
        ]));
        assert_eq!(frame.faces_half_edges, Some(ti_vec![
            vec![H(0), H(2), H(7), H(5)],
            vec![H(0), H(2), H(4), H(6)],
        ]));
        frame.assert_topology_consistent();

        // A length-2 face exists and will get squished
        let mut frame = Frame {
            ..Default::default()
        }.with_topology_vh_fh(Some(ti_vec![
            vec![H(0), H(7), H(8)],
            vec![H(1), H(2)],
            vec![H(3), H(4)],
            vec![H(5), H(6)],
            vec![H(9)],
        ]), Some(ti_vec![
            vec![H(0), H(2), H(4), H(6)],
            vec![H(8), H(9)],
        ]));
        frame.merge_two_vertices(V(0), V(4));
        assert_eq!(frame.edges_vertices, Some(ti_vec![
            [V(0), V(1)],
            [V(1), V(2)],
            [V(2), V(3)],
            [V(3), V(0)],
        ]));
        assert_eq!(frame.faces_half_edges, Some(ti_vec![
            vec![H(0), H(2), H(4), H(6)],
        ]));

        // You know what's cooler than a square cross? *2* square crosses.
        // Unfortunately, it stops being a manifold.
        let mut frame = Frame {
            frame_attributes: indexset![FrameAttribute::Manifold, FrameAttribute::Orientable],
            ..Default::default()
        }.with_topology_vh_fh(Some(ti_vec![
            vec![H(0), H(8) , H(7)],
            vec![H(2), H(10), H(1)],
            vec![H(4), H(12), H(3)],
            vec![H(6), H(14), H(5)],
            vec![H(9), H(11), H(13), H(15)],
            vec![H(16), H(24), H(23)],
            vec![H(18), H(26), H(17)],
            vec![H(20), H(28), H(19)],
            vec![H(22), H(30), H(21)],
            vec![H(25), H(27), H(29), H(31)],
        ]), Some(ti_vec![
            vec![H(0), H(10), H(9)],
            vec![H(2), H(12), H(11)],
            vec![H(4), H(14), H(13)],
            vec![H(6), H(8), H(15)],
            vec![H(16), H(26), H(25)],
            vec![H(18), H(28), H(27)],
            vec![H(20), H(30), H(29)],
            vec![H(22), H(24), H(31)],
        ]));
        frame.merge_two_vertices(V(9), V(4));
        assert_eq!(frame.frame_attributes, indexset![]);
        assert_eq!(frame.edges_vertices, Some(ti_vec![
            [V(0), V(1)],
            [V(1), V(2)],
            [V(2), V(3)],
            [V(3), V(0)],
            [V(0), V(4)],
            [V(1), V(4)],
            [V(2), V(4)],
            [V(3), V(4)],
            [V(5), V(6)],
            [V(6), V(7)],
            [V(7), V(8)],
            [V(8), V(5)],
            [V(5), V(4)],
            [V(6), V(4)],
            [V(7), V(4)],
            [V(8), V(4)],
        ]));
        assert_eq!(frame.faces_half_edges, Some(ti_vec![
            vec![H(0), H(10), H(9)],
            vec![H(2), H(12), H(11)],
            vec![H(4), H(14), H(13)],
            vec![H(6), H(8), H(15)],
            vec![H(16), H(26), H(25)],
            vec![H(18), H(28), H(27)],
            vec![H(20), H(30), H(29)],
            vec![H(22), H(24), H(31)],
        ]));
    }

    #[test]
    fn test_point_merger() {
        let positions = ti_vec![
            vector![0.0, 0.0],
            vector![0.0, 0.0],
            vector![0.0, 0.0625],
            vector![0.0, 0.125],
        ];
        // One point
        let map = PointMerger::<f64, U2>::new(
            vec![V(2)], |v| positions[v].as_view(), 0.0)
            .into_merge_map();
        assert_eq!(map, indexmap! {
            V(2) => V(2)
        });

        // Multiple points, no epsilon
        let map = PointMerger::<f64, U2>::new(
            vec![V(0), V(1), V(2), V(3)], |v| positions[v].as_view(), 0.0)
            .into_merge_map();
        assert_eq!(map, indexmap! {
            V(0) => V(0),
            V(1) => V(0),
            V(2) => V(2),
            V(3) => V(3),
        });

        // Multiple points, some epsilon
        let map = PointMerger::<f64, U2>::new(
            vec![V(0), V(1), V(2), V(3)], |v| positions[v].as_view(), 0.07)
            .into_merge_map();
        assert_eq!(map, indexmap! {
            V(0) => V(0),
            V(1) => V(0),
            V(2) => V(0),
            V(3) => V(3),
        });

        let positions = ti_vec![
            vector![0.0, 0.0],
            vector![1.0, 0.0],
            vector![0.6, 0.6],
            vector![0.0, 1.0],
            vector![-0.6, 0.6],
            vector![-1.0, 0.0],
            vector![-0.6, -0.6],
            vector![0.0, -1.0],
            vector![0.6, -0.6],
            vector![1.0, 1.0],
            vector![-1.0, 1.0],
            vector![-1.0, -1.0],
            vector![1.0, -1.0],
        ];
        // Check the neighbors
        let map = PointMerger::<f64, U2>::new(
            (0..positions.len()).map(V),
            |v| positions[v].as_view(), 1.0)
            .into_merge_map();
        assert_eq!(map, indexmap! {
            V(0) => V(0),
            V(1) => V(0),
            V(2) => V(0),
            V(3) => V(0),
            V(4) => V(0),
            V(5) => V(0),
            V(6) => V(0),
            V(7) => V(0),
            V(8) => V(0),
            V(9) => V(9),
            V(10) => V(10),
            V(11) => V(11),
            V(12) => V(12),
        });
    }
}