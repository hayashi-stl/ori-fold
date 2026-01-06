#![cfg(test)]
use std::fmt::Debug;

use exact_number::{BasedExpr, based_expr};
use nalgebra::dmatrix;
use petgraph::Graph;
use typed_index_collections::{TiSlice, TiVec, ti_vec};

use crate::filter::Coordinate;
use crate::{EdgeDatas, EdgeField, Frame, VertexDatas, VertexField, test_utils, vertices_coords};
use crate::{Vertex as V, HalfEdge as H, Face as F};

trait ToLabelVector {
    type Output: Clone + Debug + Eq;
    fn to_label_vector(&self, denominator: u32) -> Self::Output;
}

impl ToLabelVector for f64 {
    type Output = i64;
    fn to_label_vector(&self, denominator: u32) -> Self::Output {
        (self / (denominator as f64)).round() as i64
    }
}

impl ToLabelVector for BasedExpr {
    type Output = BasedExpr;
    fn to_label_vector(&self, denominator: u32) -> Self::Output {
        self.clone()
    }
}

#[derive(Clone, Debug, Eq, PartialEq, Hash, Ord, PartialOrd)]
enum PlanarWithFacesNodeLabel<T: ToLabelVector> {
    Vertex(Vec<T::Output>),
    Edge,
    Face,
}

fn planar_with_faces_graph<T: Coordinate + ToLabelVector>(frame: &Frame, denominator: u32)
    -> Graph<PlanarWithFacesNodeLabel<T>, ()>
{
    use PlanarWithFacesNodeLabel::*;

    let mut graph = Graph::new();
    let vc = vertices_coords!(<T> frame).as_ref().unwrap();
    let vh = frame.vertices_half_edges.as_ref().unwrap();
    let fh = frame.faces_half_edges.as_ref().unwrap();

    let mut vertex_nodes = vc.column_iter()
        .map(|col| graph.add_node(Vertex(col.iter().map(|c| c.to_label_vector(denominator)).collect())))
        .collect::<TiVec<V, _>>();

    let num_edges = vh.iter().flatten().copied().max().map(|h| h.0 + 1).unwrap_or(0) / 2;
    let mut edge_nodes = (0..num_edges).map(|h| graph.add_node(Edge)).collect::<TiVec<_, _>>();

    let mut face_nodes = fh.iter_enumerated()
        .map(|(f, hs)| graph.add_node(Face))
        .collect::<TiVec<_, _>>();

    // Connect vertices
    for (v, half_edges) in vh.iter_enumerated() {
        for (i, &h) in half_edges.iter().enumerate() {
            // to edge
            graph.add_edge(vertex_nodes[v], edge_nodes[h.edge()], ());
        }
    }

    // Connect faces
    for (f, half_edges) in fh.iter_enumerated() {
        for (i, &h) in half_edges.iter().enumerate() {
            // to edge
            graph.add_edge(face_nodes[f], edge_nodes[h.edge()], ());
        }
    }

    graph
}

fn assert_planar_with_faces_vh_fh_isomorphic<T: Coordinate + ToLabelVector>(
    frame: &Frame,
    expected: &Frame,
    denominator: u32,
) {
    frame.assert_topology_consistent();

    let iso = vf2::isomorphisms(
        &planar_with_faces_graph::<T>(frame, denominator),
        &planar_with_faces_graph::<T>(expected, denominator)
    )
        .default_eq()
        .first();    

    assert!(iso.is_some(), "{frame:?}\nis not isomorphic to\n{expected:?}");
}

fn planar_with_faces_test(mut frame: Frame, mut expected: Frame, epsilon: f64) {
    let denominator = (0.4 / epsilon) as u32;
    let result = frame.clone().try_into_planar_with_faces(&BasedExpr::BASELESS_ZERO)
        .unwrap();
    assert_planar_with_faces_vh_fh_isomorphic::<BasedExpr>(&result, &expected, denominator);

    // and now to fuzz
    frame.calc_approx_coordinates();
    expected.calc_approx_coordinates();
    let coords = frame.vertices_coords_f64.as_ref().unwrap().clone();
    let mut rng = test_utils::new_rng();
    for i in 0..10 {
        println!("Iteration {i} (if 0, no fuzz)");
        let result = frame.clone().try_into_planar_with_faces(&epsilon)
            .unwrap_or_else(|e| panic!("failed at iter. {i} with {e:?}"));
        assert_planar_with_faces_vh_fh_isomorphic::<f64>(&result, &expected, denominator);
        frame.vertices_coords_f64 = Some(test_utils::fuzz_coordinates(coords.clone(), epsilon / 2.0, &mut rng))
    }
}

#[test]
fn test_try_into_planar_with_faces_line() {
    // A line to check the sanity
    let mut frame = Frame::new().init_topology_coords_exact_ev(dmatrix![
        based_expr!(0), based_expr!(0);
        based_expr!(1), based_expr!(0);
    ].transpose(), ti_vec![
        [V(0), V(1)],
    ]);
    let mut expected = frame.clone();
    expected.init_face_data([]);
    planar_with_faces_test(frame.clone(), expected.clone(), 1e-8 );
}

#[test]
fn test_try_into_planar_with_faces_cross() {
    // Now an intersection occurs
    let frame = Frame::new().init_topology_coords_exact_ev(dmatrix![
        based_expr!(0), based_expr!(0);
        based_expr!(1), based_expr!(0);
        based_expr!(1/2), based_expr!(1/2);
        based_expr!(1/2), based_expr!(-1/2);
    ].transpose(), ti_vec![
        [V(0), V(1)],
        [V(2), V(3)],
    ]);
    let mut expected = Frame::new().init_topology_coords_exact_ev(dmatrix![
        based_expr!(0), based_expr!(0);
        based_expr!(1), based_expr!(0);
        based_expr!(1/2), based_expr!(1/2);
        based_expr!(1/2), based_expr!(-1/2);
        based_expr!(1/2), based_expr!(0);
    ].transpose(), ti_vec![
        [V(0), V(4)],
        [V(1), V(4)],
        [V(2), V(4)],
        [V(3), V(4)],
    ]).with_topology_fh(ti_vec![]);
    planar_with_faces_test(frame.clone(), expected.clone(), 1e-8 );
}

#[test]
fn test_try_into_planar_with_faces_square() {
    let frame = Frame::new().init_topology_coords_exact_ev(dmatrix![
        based_expr!(0), based_expr!(0);
        based_expr!(1), based_expr!(0);
        based_expr!(0), based_expr!(1);
        based_expr!(1), based_expr!(1);
    ].transpose(), ti_vec![
        [V(0), V(1)],
        [V(1), V(3)],
        [V(3), V(2)],
        [V(2), V(0)],
    ]);
    let expected = Frame::new().init_topology_coords_exact_ev(dmatrix![
        based_expr!(0), based_expr!(0);
        based_expr!(1), based_expr!(0);
        based_expr!(0), based_expr!(1);
        based_expr!(1), based_expr!(1);
    ].transpose(), ti_vec![
        [V(0), V(1)],
        [V(1), V(3)],
        [V(3), V(2)],
        [V(2), V(0)],
    ]).with_topology_fh(ti_vec![
        vec![H(0), H(2), H(4), H(6)]
    ]);
    planar_with_faces_test(frame.clone(), expected.clone(), 1e-8 );
}

#[test]
fn test_try_into_planar_with_faces_broken_square() {
    let frame = Frame::new().init_topology_coords_exact_ev(dmatrix![
        based_expr!(0), based_expr!(0);
        based_expr!(0), based_expr!(0);
        based_expr!(1), based_expr!(0);
        based_expr!(1), based_expr!(0);
        based_expr!(0), based_expr!(1);
        based_expr!(0), based_expr!(1);
        based_expr!(1), based_expr!(1);
        based_expr!(1), based_expr!(1);
    ].transpose(), ti_vec![
        [V(0), V(2)],
        [V(3), V(6)],
        [V(7), V(4)],
        [V(5), V(1)],
    ]);
    let expected = Frame::new().init_topology_coords_exact_ev(dmatrix![
        based_expr!(0), based_expr!(0);
        based_expr!(1), based_expr!(0);
        based_expr!(0), based_expr!(1);
        based_expr!(1), based_expr!(1);
    ].transpose(), ti_vec![
        [V(0), V(1)],
        [V(1), V(3)],
        [V(3), V(2)],
        [V(2), V(0)],
    ]).with_topology_fh(ti_vec![
        vec![H(0), H(2), H(4), H(6)]
    ]);
    planar_with_faces_test(frame.clone(), expected.clone(), 1e-8 );
}

#[test]
fn test_try_into_planar_with_faces_squares() {
    // +-----+
    // |     |
    // |  +--+--+
    // |  |  |  |
    // +--+--+  |
    //    |     |
    //    +-----+
    let frame = Frame::new().init_topology_coords_exact_ev(dmatrix![
        based_expr!(0), based_expr!(1);
        based_expr!(2), based_expr!(1);
        based_expr!(0), based_expr!(3);
        based_expr!(2), based_expr!(3);
        based_expr!(1), based_expr!(0);
        based_expr!(3), based_expr!(0);
        based_expr!(1), based_expr!(2);
        based_expr!(3), based_expr!(2);
    ].transpose(), ti_vec![
        [V(0), V(1)],
        [V(1), V(3)],
        [V(3), V(2)],
        [V(2), V(0)],
        [V(4), V(5)],
        [V(5), V(7)],
        [V(7), V(6)],
        [V(6), V(4)],
    ]);
    let expected = Frame::new().init_topology_coords_exact_ev(dmatrix![
        based_expr!(0), based_expr!(1);
        based_expr!(2), based_expr!(1);
        based_expr!(0), based_expr!(3);
        based_expr!(2), based_expr!(3);
        based_expr!(1), based_expr!(0);
        based_expr!(3), based_expr!(0);
        based_expr!(1), based_expr!(2);
        based_expr!(3), based_expr!(2);
        based_expr!(1), based_expr!(1);
        based_expr!(2), based_expr!(2);
    ].transpose(), ti_vec![
        [V(0), V(8)],
        [V(8), V(1)],
        [V(1), V(9)],
        [V(9), V(3)],
        [V(3), V(2)],
        [V(2), V(0)],
        [V(4), V(5)],
        [V(5), V(7)],
        [V(7), V(9)],
        [V(9), V(6)],
        [V(6), V(8)],
        [V(8), V(4)],
    ]).with_topology_fh(ti_vec![
        vec![H(0), H(21), H(19), H(6), H(8), H(10)],
        vec![H(16), H(5), H(3), H(22), H(12), H(14)],
        vec![H(2), H(4), H(18), H(20)],
    ]);
    planar_with_faces_test(frame.clone(), expected.clone(), 1e-8 );
}
