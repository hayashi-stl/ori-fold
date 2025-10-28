use std::mem;

use approx::{relative_eq, relative_ne};
use indexmap::IndexMap;
use nalgebra::{vector, Affine2, ClosedSubAssign, DMatrix, DMatrixView, Matrix2xX, Matrix3, Point2, RealField, Reflection2, Scalar, Vector2};
use num_traits::RefNum;
use typed_index_collections::{ti_vec, TiSlice, TiVec};

use crate::{filter, fold::{EdgesVerticesEx, CoordsRef, Edge, EdgeAssignment, EdgesFaceCornersEx, EdgesFaceCornersSlice, EdgesVerticesSlice, Face, FaceCorner, FacesHalfEdgesSlice, Fold, Frame, FrameAttribute, HalfEdge, Vertex}, geom::{sort_by_angle_ref, Atan2, MatrixView2Dyn, NumEx}, manifold::OrientableError};
use crate::geom;

/// Assuming a locally flat foldable crease pattern in the xy plane
/// (and thus an orientable manifold) without length-0 edges,
/// computes:
/// * the flat-folded geometry as determined by repeated reflection relative to `root_face`,
/// * the transformation matrix mapping each face's unfolded --> folded geometry
/// * for each face, `true` if the face got reflected and `false` otherwise.
/// 
/// For connected-components-by-face outside of the one in `root_face`:
/// * the flat-folded geometry is the original geometry
/// * the face transforms are the identity matrix
/// * the face doesn't get reflected, so `false` for orientation
pub fn flat_folded_geometry<T: NumEx + RealField>(
    vertices_coords: MatrixView2Dyn<T>,
    edges_vertices: &EdgesVerticesSlice,
    edges_face_corners: &EdgesFaceCornersSlice,
    faces_half_edges: &FacesHalfEdgesSlice,
    root_face: Face,
) -> (Matrix2xX<T>, TiVec<Face, Affine2<T>>, TiVec<Face, bool>)
    where for<'a> &'a T: RefNum<T>
{
    //let mut max_error2 = T::zero();
    let mut level = vec![root_face];

    let mut faces_flat_fold_transform = ti_vec![Affine2::from_matrix_unchecked(Matrix3::identity()); faces_half_edges.len()];
    let mut faces_flat_fold_orientation = ti_vec![false; faces_half_edges.len()];
    let mut vertices_flat_fold_coords = vertices_coords.clone().into_owned();
    let mut faces_reached = ti_vec![false; faces_half_edges.len()];
    let mut vertices_reached = ti_vec![false; vertices_coords.ncols()];
    faces_reached[root_face] = true;

    while !level.is_empty() {
        let mut next_level = vec![];

        for face in level {
            let orientation = !faces_flat_fold_orientation[face];
            for &half_edge in &faces_half_edges[face] {
                if let Some(&FaceCorner(face_2, _)) = edges_face_corners.at(half_edge.flipped()).first() {
                    let coords = edges_vertices.at(half_edge).map(|v| vertices_coords.column(v.0));
                    let reflection = geom::reflect_line_matrix(coords[0].as_view(), coords[1].as_view());
                    let transform = &faces_flat_fold_transform[face] * reflection;

                    if mem::replace(&mut faces_reached[face_2], true) {
                        //let diff = transform.matrix() - faces_flat_fold_transform[face_2].matrix();
                        //max_error2 = max_error2
                        //    .max(diff.fixed_view::<2, 1>(0, 0).norm_squared())
                        //    .max(diff.fixed_view::<2, 1>(0, 1).norm_squared())
                        //    .max(diff.fixed_view::<2, 1>(0, 2).norm_squared());
                    } else {
                        faces_flat_fold_transform[face_2] = transform;
                        faces_flat_fold_orientation[face_2] = orientation;
                        let transform = &faces_flat_fold_transform[face_2];

                        for &half_edge_2 in &faces_half_edges[face_2] {
                            for vertex_2 in edges_vertices.at(half_edge_2) {
                                let mapped = geom::transform(transform, vertices_coords.column(vertex_2.0).as_view());
                                if mem::replace(&mut vertices_reached[vertex_2], true) {
                                    //max_error2 = max_error2.max((vertices_flat_fold_coords.column(vertex_2.0) - mapped).norm_squared())
                                } else {
                                    vertices_flat_fold_coords.set_column(vertex_2.0, &mapped);
                                }
                            }
                        }
                        next_level.push(face_2);
                    }
                }
            }
        }
        level = next_level;
    }
    
    (vertices_flat_fold_coords, faces_flat_fold_transform, faces_flat_fold_orientation)
}

/// An error resulting from trying to flat-fold geometry without taking layering into account
#[derive(Clone, Debug)]
#[cfg_attr(test, derive(PartialEq, Eq, PartialOrd, Ord))] // Really no point in this derivation outside of tests
pub enum FlatFoldError {
    NotOrientable(OrientableError),
    Not2D { dims: usize },
    NoDefinedFaces,
    ZeroLengthEdge { edge: Edge, vertices: [Vertex; 2] },
    VertexNotClosed { vertex: Vertex },
}

impl Frame {
    /// Computes the flat-folded geometry as in `Frame::try_flat_folded_geometry` without
    /// doing any checks. Beware.
    /// 
    /// In particular, the following must be true:
    /// * the frame is not just orientable, but *oriented*, with `FrameAttribute::Orientable` set
    /// * the coordinates are 2D
    /// * `faces_half_edges` exists (this is at least temporary)
    /// * no edges have length 0
    /// * all vertices are closed according to Kawasaki's theorem
    pub fn flat_folded_geometry_unchecked<T: NumEx + RealField>(&self, root_face: Face, vertices_coords: DMatrixView<'_, T>)
        -> (Matrix2xX<T>, TiVec<Face, Affine2<T>>, TiVec<Face, bool>) where
        for<'a> &'a T: RefNum<T>
    {
        let vertices_coords = vertices_coords.fixed_rows::<2>(0);
        let faces_half_edges = self.faces_half_edges.as_ref().unwrap();
        let edges_face_corners = self.edges_face_corners.as_ref().unwrap();
        let edges_vertices = self.edges_vertices.as_ref().unwrap();
        flat_folded_geometry(vertices_coords, edges_vertices, edges_face_corners, faces_half_edges, root_face)
    }

    /// Computes:
    /// * the flat-folded geometry as determined by repeated reflection relative to `root_face`,
    /// * the transformation matrix mapping each face's unfolded --> folded geometry
    /// * for each face, `true` if the face got reflected and `false` otherwise.
    /// 
    /// For connected-components-by-face outside of the one in `root_face`:
    /// * the flat-folded geometry is the original geometry
    /// * the face transforms are the identity matrix
    /// * the face doesn't get reflected, so `false` for orientation
    /// 
    /// # Errors
    /// Returns an error if
    /// * the frame is not orientable
    /// * the coordinates are not 2D
    /// * `faces_half_edges` does not exist (this is at least temporary)
    /// * an edge has length 0
    /// * a vertex isn't closed according to Kawasaki's theorem
    pub fn try_flat_folded_geometry<T: NumEx + RealField>(&self, root_face: Face, vertices_coords: DMatrixView<'_, T>)
        -> Result<(Matrix2xX<T>, TiVec<Face, Affine2<T>>, TiVec<Face, bool>), Vec<FlatFoldError>> where
        for<'a> &'a T: RefNum<T>
    {
        let orientable = if self.frame_attributes.contains(&FrameAttribute::Orientable) {
            self
        } else {
            &self.clone().try_into_orientable()
                .map_err(|e| e.into_iter().map(FlatFoldError::NotOrientable).collect::<Vec<_>>())?.0
        };

        if vertices_coords.nrows() != 2 {
            Err(vec![FlatFoldError::Not2D { dims: vertices_coords.nrows() }])?
        }
        let vertices_coords = vertices_coords.fixed_rows::<2>(0);
        let faces_half_edges = orientable.faces_half_edges.as_ref().ok_or(vec![FlatFoldError::NoDefinedFaces])?;
        let edges_face_corners = orientable.edges_face_corners.as_ref().unwrap(); // guaranteed to exist if the above exists
        let edges_vertices = orientable.edges_vertices.as_ref().unwrap(); // guaranteed to exist if the above exists
        let vertices_half_edges = orientable.vertices_half_edges.as_ref().unwrap(); // guaranteed to exist if the above exists
        let mut errors = vec![];

        // Edge lengths
        for (e, &vertices) in edges_vertices.iter_enumerated() {
            let diff = vertices_coords.column(vertices[0].0) - vertices_coords.column(vertices[1].0);
            if relative_eq!(diff.norm_squared(), T::zero()) {
                errors.push(FlatFoldError::ZeroLengthEdge { edge: e, vertices: vertices })
            }
        }
        if !errors.is_empty() { return Err(errors); }

        // Check closure
        for (v, half_edges) in vertices_half_edges.iter_enumerated() {
            // Only loops matter
            if half_edges.iter().any(|&h| edges_face_corners.at(h).is_empty()) { continue }

            let mut even_faces = true;
            let mut vector = vector![T::one(), T::zero()];
            let coords = vertices_coords.column(v.0);
            for &h in half_edges {
                let other = vertices_coords.column(edges_vertices.at(h)[1].0);
                vector = geom::reflect(vector.as_view(), (other - &coords).as_view()); // already checked for 0-length vectors
                even_faces = !even_faces;
            }

            if !even_faces || relative_ne!(vector, vector![T::one(), T::zero()]) {
                errors.push(FlatFoldError::VertexNotClosed { vertex: v })
            }
        }
        if !errors.is_empty() { return Err(errors); }

        Ok(flat_folded_geometry(vertices_coords, edges_vertices, edges_face_corners, faces_half_edges, root_face))
    }

}

#[cfg(test)]
mod test {
    use exact_number::based_expr;
    use nalgebra::{Affine2, DMatrix, DMatrixView, Matrix2xX, RealField};
    use num_traits::RefNum;
    use typed_index_collections::{ti_vec, TiSlice, TiVec};

    use crate::{fold::{Face as F, Frame, HalfEdge as H}, geom::NumEx};

    fn flat_folded_geometry_success_test<T: NumEx + RealField>(
        frame: Frame,
        root_face: F,
        coords: DMatrix<T>,
        expected_geometry: Matrix2xX<T>,
        expected_transforms: TiVec<F, Affine2<T>>,
        expected_reflecteds: TiVec<F, bool>
    ) where 
        for<'a> &'a T: RefNum<T>
    {
        let result = (expected_geometry, expected_transforms, expected_reflecteds);
        assert_eq!(frame.flat_folded_geometry_unchecked(root_face, coords.as_view()), result);
        assert_eq!(frame.try_flat_folded_geometry(root_face, coords.as_view()), Ok(result));
    }

    macro_rules! exact_affine_2 {
        (($($a:tt)*), ($($b:tt)*), ($($c:tt)*); ($($d:tt)*), ($($e:tt)*), ($($f:tt)*);) => {
            nalgebra::Affine2::from_matrix_unchecked(nalgebra::matrix![
                based_expr!($($a)*), based_expr!($($b)*), based_expr!($($c)*);
                based_expr!($($d)*), based_expr!($($e)*), based_expr!($($f)*);
                exact_number::BasedExpr::BASELESS_ZERO, exact_number::BasedExpr::BASELESS_ZERO, exact_number::BasedExpr::BASELESS_ONE;
            ])
        };
    }

    #[test]
    fn test_flat_folded_geometry() {
        // A simple square.
        flat_folded_geometry_success_test(Frame {..Default::default()}.with_topology_vh_fh(Some(ti_vec![
            vec![H(0), H(7)],
            vec![H(2), H(1)],
            vec![H(4), H(3)],
            vec![H(6), H(5)],
        ]), Some(ti_vec![
            vec![H(0), H(2), H(4), H(6)],
        ])), F(0), DMatrix::from_vec(2, 4, vec![
            based_expr!(0), based_expr!(0),
            based_expr!(1), based_expr!(0),
            based_expr!(1), based_expr!(1),
            based_expr!(0), based_expr!(1),
        ]), Matrix2xX::from_vec(vec![
            based_expr!(0), based_expr!(0),
            based_expr!(1), based_expr!(0),
            based_expr!(1), based_expr!(1),
            based_expr!(0), based_expr!(1),
        ]), ti_vec![
            exact_affine_2![(1), (0), (0); (0), (1), (0);]
        ], ti_vec![
            false
        ]);

        // A square attached to a triangle.
        flat_folded_geometry_success_test(Frame {..Default::default()}.with_topology_vh_fh(Some(ti_vec![
            vec![H(0), H(7)],
            vec![H(8), H(2), H(1)],
            vec![H(4), H(3), H(11)],
            vec![H(6), H(5)],
            vec![H(10), H(9)],
        ]), Some(ti_vec![
            vec![H(0), H(2), H(4), H(6)],
            vec![H(8), H(10), H(3)],
        ])), F(0), DMatrix::from_vec(2, 5, vec![
            based_expr!(0), based_expr!(0),
            based_expr!(1), based_expr!(0),
            based_expr!(1), based_expr!(1),
            based_expr!(0), based_expr!(1),
            based_expr!(7/4), based_expr!(1/2),
        ]), Matrix2xX::from_vec(vec![
            based_expr!(0), based_expr!(0),
            based_expr!(1), based_expr!(0),
            based_expr!(1), based_expr!(1),
            based_expr!(0), based_expr!(1),
            based_expr!(1/4), based_expr!(1/2),
        ]), ti_vec![
            exact_affine_2![(1), (0), (0); (0), (1), (0);],
            exact_affine_2![(-1), (0), (2); (0), (1), (0);],
        ], ti_vec![
            false,
            true,
        ]);
    }
}