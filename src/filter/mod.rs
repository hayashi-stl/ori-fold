use std::{iter, mem};

use exact_number::BasedExpr;
use indexmap::{indexmap, map::Entry, IndexMap};
use nalgebra::{DMatrix, DVector, Dyn, Scalar};
use typed_index_collections::{ti_vec, TiVec};

use crate::{EdgeDatas, FaceDatas, VertexDatas, filter::intersect::IntersectCoordinate, fold::{AtFaceCorner, CoordsRef, Edge, EdgeAssignment, EdgeData, EdgeField, EdgesFaceCornersEx, EdgesVerticesEx, Face, FaceCorner, FaceData, Frame, FrameAttribute, HalfEdge, Vertex, VertexData, VertexField}};

pub mod intersect;
pub mod split_merge;

#[derive(Clone, Copy)]
struct SwapRemove;
#[derive(Clone, Copy)]
struct ShiftRemove;
trait RemoveStrategy: Clone {
    fn index_map<K>(&self, vec_len: usize, key: K, k: K) -> Option<K> where
        K: PartialOrd + From<usize>,
        usize: From<K>;
    fn remove<K, T>(&self, vec: &mut TiVec<K, T>, key: K) -> T where usize: From<K>;
    fn remove_column<T: Scalar>(&self, vec: &mut DMatrix<T>, column: usize) -> DVector<T>;
}

impl RemoveStrategy for SwapRemove {
    fn index_map<K>(&self, vec_len: usize, key: K, k: K) -> Option<K> where
        K: PartialOrd + From<usize>,
        usize: From<K>
    {
        if k == key { None }
        else if k == K::from(vec_len - 1) { Some(key) }
        else { Some(k) }
    }

    fn remove<K, T>(&self, vec: &mut TiVec<K, T>, key: K) -> T where usize: From<K> { vec.swap_remove(key) }

    fn remove_column<T: Scalar>(&self, vec: &mut DMatrix<T>, column: usize) -> DVector<T> {
        let result = vec.column(column).into_owned();
        vec.swap_columns(column, vec.ncols() - 1);
        let ncols = vec.ncols();
        *vec = mem::take(vec).remove_column(ncols - 1);
        result
    }
}

impl RemoveStrategy for ShiftRemove {
    fn index_map<K>(&self, _vec_len: usize, key: K, k: K) -> Option<K> where
        K: PartialOrd + From<usize>,
        usize: From<K>
    {
        if k == key { None }
        else if k < key { Some(k) }
        else { Some(K::from(usize::from(k) - 1)) }
    }

    fn remove<K, T>(&self, vec: &mut TiVec<K, T>, key: K) -> T where usize: From<K> { vec.remove(key) }

    fn remove_column<T: Scalar>(&self, vec: &mut DMatrix<T>, column: usize) -> DVector<T> {
        let result = vec.column(column).into_owned();
        *vec = mem::take(vec).remove_column(column);
        result
    }
}

pub trait Coordinate: Sized + IntersectCoordinate {
    fn vertices_coords(frame: &'_ Frame) -> &'_ Option<DMatrix<Self>>;
    fn vertices_coords_mut(frame: &'_ mut Frame) -> &'_ mut Option<DMatrix<Self>>;
}

impl Coordinate for f64 {
    fn vertices_coords(frame: &'_ Frame) -> &'_ Option<DMatrix<Self>> { &frame.vertices_coords_f64 }
    fn vertices_coords_mut(frame: &'_ mut Frame) -> &'_ mut Option<DMatrix<Self>> { &mut frame.vertices_coords_f64 }
}

impl Coordinate for BasedExpr {
    fn vertices_coords(frame: &'_ Frame) -> &'_ Option<DMatrix<Self>> { &frame.vertices_coords_exact }
    fn vertices_coords_mut(frame: &'_ mut Frame) -> &'_ mut Option<DMatrix<Self>> { &mut frame.vertices_coords_exact }
}

/// Sets a column without cloning.
fn set_column<T: Scalar>(mtx: &mut DMatrix<T>, index: usize, col: DVector<T>) {
    let dims = col.nrows();
    if col.nrows() != mtx.nrows() {
        panic!("dimension mismatch: tried to set a {}D entry in a {}D vector array", col.nrows(), mtx.nrows())
    }
    let slice = mtx.column_mut(index).data.into_slice_mut();
    let mut col: Vec<T> = col.data.into();
    slice.swap_with_slice(&mut col);
}

fn insert_column<T: Scalar>(mtx: &mut DMatrix<T>, col: DVector<T>) {
    let dims = col.nrows();
    insert_columns(mtx, col.reshape_generic(Dyn(dims), Dyn(1)))
}

fn insert_columns<T: Scalar>(mtx: &mut DMatrix<T>, cols: DMatrix<T>) {
    let index = mtx.ncols();
    let len = cols.ncols();
    let dims = cols.nrows();
    if cols.nrows() != mtx.nrows() {
        panic!("dimension mismatch: tried to insert {}D entries to a {}D vector array", cols.nrows(), mtx.nrows())
    }
    unsafe {
        let mut uninit = mem::take(mtx).insert_columns_generic_uninitialized(index, Dyn(len));
        let mut cols = cols.data.resize(len * dims).into_iter();
        for c in 0..len {
            for r in 0..dims {
                uninit[(r, index + c)].write(cols.next().unwrap().assume_init());
            }
        }
        // Safety: All new entries are now initialized
        *mtx = uninit.assume_init();
    }
}

impl Frame {
    /// Gets a reference to the vertex coordinates, either exact or approximate.
    pub fn coords_ref(&self) -> Option<CoordsRef<'_>> {
        if let Some(coords) = self.vertices_coords_exact.as_ref() {
            Some(CoordsRef::Exact(coords))
        } else if let Some(coords) = self.vertices_coords_f64.as_ref() {
            Some(CoordsRef::Approx(coords))
        } else { None }
    }

    fn set_num_dimensions(&mut self, num_dimensions: usize) {
        self.frame_attributes.swap_remove(&FrameAttribute::Abstract);
        self.frame_attributes.swap_remove(&FrameAttribute::_2D);
        self.frame_attributes.swap_remove(&FrameAttribute::_3D);
        match num_dimensions {
            2 => self.frame_attributes.insert(FrameAttribute::_2D),
            3 => self.frame_attributes.insert(FrameAttribute::_3D),
            _ => true
        };
    }

    /// Initializes the vertex data with some fields. This should be called only on a frame with no vertex data,
    /// and before adding new vertices.
    /// 
    /// `vertices_half_edges` does not get initialized.
    /// 
    /// The coordinates will have `num_dimensions` dimensions if they are specified to exist.
    pub fn init_vertex_data(&mut self, fields: impl IntoIterator<Item = VertexField>, num_dimensions: usize) {
        assert_eq!(self.num_vertices, 0);
        assert_eq!(self.vertices_coords_f64, None);
        assert_eq!(self.vertices_coords_exact, None);
        assert_eq!(self.vertices_half_edges, None);
        assert_eq!(self.vertices_custom.len(), 0);

        for field in fields.into_iter() {
            match field {
                VertexField::CoordsF64 => {
                    self.set_num_dimensions(num_dimensions);
                    self.vertices_coords_f64 = Some(DMatrix::from_vec(num_dimensions, 0, vec![]))
                },
                VertexField::CoordsExact => {
                    self.set_num_dimensions(num_dimensions);
                    self.vertices_coords_exact = Some(DMatrix::from_vec(num_dimensions, 0, vec![]))
                },
                VertexField::HalfEdges => (),
                VertexField::Custom(k) => { self.vertices_custom.insert(k, ti_vec![]); }
            }
        }
    }

    /// Initializes the edge data with some fields. This should be called only on a frame with no edge data,
    /// and before adding new edges.
    /// 
    /// `vertices_half_edges` and `edges_vertices` always get initialized.
    /// `edges_face_corners` does not get initialized.
    pub fn init_edge_data(&mut self, fields: impl IntoIterator<Item = EdgeField>) {
        assert_eq!(self.edges_vertices, None);
        assert_eq!(self.edges_face_corners, None);
        assert_eq!(self.edges_assignment, None);
        assert_eq!(self.edges_fold_angle_f64, None);
        assert_eq!(self.edges_fold_angle_exact, None);
        assert_eq!(self.edges_length_f64, None);
        assert_eq!(self.edges_length2_exact, None);
        assert_eq!(self.edges_custom.len(), 0);

        self.vertices_half_edges = Some(ti_vec![]);
        self.edges_vertices = Some(ti_vec![]);
        for field in fields.into_iter() {
            match field {
                EdgeField::Vertices => (),
                EdgeField::FaceCorners => (),
                EdgeField::Assignment => self.edges_assignment = Some(ti_vec![]),
                EdgeField::FoldAngleF64 => self.edges_fold_angle_f64 = Some(ti_vec![]),
                EdgeField::FoldAngleExact => self.edges_fold_angle_exact = Some(ti_vec![]),
                EdgeField::LengthF64 => self.edges_length_f64 = Some(ti_vec![]),
                EdgeField::Length2Exact => self.edges_length2_exact = Some(ti_vec![]),
                EdgeField::Custom(k) => { self.edges_custom.insert(k, ti_vec![]); }
            };
        }
    }

    /// Adds a vertex given some vertex data and returns its index.
    /// 
    /// All fields that exists in this frame must have their respective field specified
    /// in `vertex_data`
    /// 
    /// Note that `vertex_data.half_edges` gets ignored.
    pub fn add_vertex(&mut self, mut vertex_data: VertexData) -> Vertex {
        let index = Vertex(self.num_vertices);
        self.num_vertices += 1;

        self.vertices_coords_f64.as_mut().map(|vec| {
            insert_column(vec,  vertex_data.coords_f64.take().unwrap());
        });
        self.vertices_coords_exact.as_mut().map(|vec| {
            insert_column(vec,  vertex_data.coords_exact.take().unwrap());
        });
        self.vertices_half_edges.as_mut().map(|vec| vec.push(vec![]));

        self.vertices_custom.iter_mut().for_each(|(k, vec)| vec.push(vertex_data.custom.swap_remove(k).unwrap()));

        index
    }

    /// Adds vertices given some vertex data and returns the index of the first new vertex.
    /// 
    /// All fields that exists in this frame must have their respective field specified
    /// in `vertex_data`, and all fields must have the same number of elements.
    /// 
    /// Note that `vertex_data.half_edges` gets ignored.
    pub fn add_vertices(&mut self, mut vertex_datas: VertexDatas) -> Vertex {
        let index = Vertex(self.num_vertices);
        let len = vertex_datas.len();
        self.num_vertices += len;

        self.vertices_coords_f64.as_mut().map(|vec| {
            insert_columns(vec, vertex_datas.coords_f64.take().unwrap());
        });
        self.vertices_coords_exact.as_mut().map(|vec| {
            insert_columns(vec, vertex_datas.coords_exact.take().unwrap());
        });
        self.vertices_half_edges.as_mut().map(|vec| vec.push(vec![]));

        self.vertices_custom.iter_mut().for_each(|(k, vec)| vec.extend(vertex_datas.custom.swap_remove(k).unwrap()));

        index
    }

    /// Gets the data for a vertex all in one place
    pub fn vertex_data(&mut self, vertex: Vertex) -> VertexData {
        VertexData {
            coords_f64: self.vertices_coords_f64.as_ref().map(|vec| vec.column(vertex.0).into_owned()),
            coords_exact: self.vertices_coords_exact.as_ref().map(|vec| vec.column(vertex.0).into_owned()),
            half_edges: self.vertices_half_edges.as_ref().map(|vec| vec[vertex].clone()),
            custom: self.vertices_custom.iter_mut().map(|(k, vec)| (k.clone(), vec[vertex].clone())).collect()
        }
    }

    /// Adds a vertex with the same data as the vertex given,
    /// except at position `at`.
    /// If vertex coordinates with type `T` do not exist, `at` gets ignored.
    /// 
    /// Note that `vertex_data.half_edges` gets ignored.
    pub fn add_vertex_like<T: Coordinate>(&mut self, like: Vertex, at: DVector<T>) -> Vertex {
        let data = self.vertex_data(like);
        let v = self.add_vertex(data);
        T::vertices_coords_mut(self).as_mut().map(|vc| set_column(vc, v.0, at));
        v
    }

    /// Adds an edge given some edge data and returns its index
    /// without checking whether it preserves the frame attributes.
    /// If `fix_vertices_half_edges` is false, then `vertices_half_edges` is left unmodified.
    /// See [`add_edge`](Frame::add_edge) for further documentation.
    pub fn add_edge_unchecked(&mut self, mut edge_data: EdgeData, fix_vertices_half_edges: bool) -> Edge {
        let index = Edge(self.edges_vertices.as_ref().unwrap().len());

        self.edges_vertices.as_mut().unwrap().push(edge_data.vertices);
        if fix_vertices_half_edges {
            self.vertices_half_edges.as_mut().unwrap()[edge_data.vertices[0]].push(HalfEdge::new(index, false));
            self.vertices_half_edges.as_mut().unwrap()[edge_data.vertices[1]].push(HalfEdge::new(index, true));
        }
        self.edges_face_corners.as_mut().map(|vec| vec.push([vec![], vec![]]));
        self.edges_assignment.as_mut().map(|vec| vec.push(edge_data.assignment.unwrap()));
        self.edges_fold_angle_f64.as_mut().map(|vec| vec.push(edge_data.fold_angle_f64.unwrap()));
        self.edges_fold_angle_exact.as_mut().map(|vec| vec.push(edge_data.fold_angle_exact.unwrap()));
        self.edges_length_f64.as_mut().map(|vec| vec.push(edge_data.length_f64.unwrap()));
        self.edges_length2_exact.as_mut().map(|vec| vec.push(edge_data.length2_exact.unwrap()));
        self.edges_custom.iter_mut().for_each(|(k, vec)| vec.push(edge_data.custom.swap_remove(k).unwrap()));
        index
    }

    /// Adds edges given some edge datas and returns the index of the first edge
    /// without checking whether adding the edges preserves the frame attributes.
    /// See [`add_edges`](Frame::add_edges) for further documentation.
    pub fn add_edges_unchecked(&mut self, mut edge_datas: EdgeDatas) -> Edge {
        let index = Edge(self.edges_vertices.as_ref().unwrap().len());
        let len = edge_datas.len();

        for (i, &[v0, v1]) in edge_datas.vertices.iter().enumerate() {
            self.vertices_half_edges.as_mut().unwrap()[v0].push(HalfEdge::new(index + Edge(i), false));
            self.vertices_half_edges.as_mut().unwrap()[v1].push(HalfEdge::new(index + Edge(i), true));
        }
        self.edges_vertices.as_mut().unwrap().extend(edge_datas.vertices);
        self.edges_face_corners.as_mut().map(|vec| vec.extend(vec![[vec![], vec![]]; len]));
        self.edges_assignment.as_mut().map(|vec| vec.extend(edge_datas.assignment.unwrap()));
        self.edges_fold_angle_f64.as_mut().map(|vec| vec.extend(edge_datas.fold_angle_f64.unwrap()));
        self.edges_fold_angle_exact.as_mut().map(|vec| vec.extend(edge_datas.fold_angle_exact.unwrap()));
        self.edges_length_f64.as_mut().map(|vec| vec.extend(edge_datas.length_f64.unwrap()));
        self.edges_length2_exact.as_mut().map(|vec| vec.extend(edge_datas.length2_exact.unwrap()));
        self.edges_custom.iter_mut().for_each(|(k, vec)| vec.extend(edge_datas.custom.swap_remove(k).unwrap()));
        index
    }

    /// Adds an edge given some edge data and returns its index.
    /// 
    /// * If the edge is a [`Cut`](EdgeAssignment::Cut) then [`NoCuts`](FrameAttribute::NoCuts) is cleared.
    /// * If the edge is a [`Join`](EdgeAssignment::Join) then [`NoJoins`](FrameAttribute::NoJoins) is cleared.
    /// * If the mesh is no longer a manifold, then [`Manifold`](FrameAttribute::Manifold) is cleared.
    /// 
    /// All fields that exists in this frame must have their respective field specified
    /// in `edge_data`
    /// 
    /// Note that `edge_data.face_corners` gets ignored.
    pub fn add_edge(&mut self, edge_data: EdgeData) -> Edge {
        let [v0, v1] = edge_data.vertices;
        assert_ne!(v0, v1, "self-loops are not allowed, not even between {v0} and {v1}");

        if self.frame_attributes.contains(&FrameAttribute::NoCuts) && edge_data.assignment == Some(EdgeAssignment::Cut) {
            self.frame_attributes.swap_remove(&FrameAttribute::NoCuts);
        }
        if self.frame_attributes.contains(&FrameAttribute::NoJoins) && edge_data.assignment == Some(EdgeAssignment::Join) {
            self.frame_attributes.swap_remove(&FrameAttribute::NoJoins);
        }
        if self.frame_attributes.contains(&FrameAttribute::Manifold) {
            let vh = self.vertices_half_edges.as_ref().unwrap();
            let ec = self.edges_face_corners.as_ref().unwrap();
            for v in edge_data.vertices {
                if vh[v].iter().all(|&h| ec.at(h).len() + ec.at(h.flipped()).len() == 2) {
                    self.frame_attributes.swap_remove(&FrameAttribute::Orientable);
                    self.frame_attributes.swap_remove(&FrameAttribute::Manifold);
                    break;
                }
            }
        }
        self.add_edge_unchecked(edge_data, true)
    }

    /// Adds multiple edges given some edge data and returns the index of the first edge.
    /// 
    /// If the mesh is no longer a manifold, then [`Manifold`](FrameAttribute::Manifold) is cleared.
    /// 
    /// All fields that exists in this frame must have their respective field specified
    /// in `edge_data`
    /// 
    /// Note that `edge_data.face_corners` gets ignored.
    pub fn add_edges(&mut self, edge_datas: EdgeDatas) -> Edge {
        for &[v0, v1] in &edge_datas.vertices {
            assert_ne!(v0, v1, "self-loops are not allowed, not even between {v0} and {v1}");
        }

        if self.frame_attributes.contains(&FrameAttribute::NoCuts) &&
            edge_datas.assignment.as_ref().map(|vec| vec.contains(&EdgeAssignment::Cut)).unwrap_or(false)
        {
            self.frame_attributes.swap_remove(&FrameAttribute::NoCuts);
        }
        if self.frame_attributes.contains(&FrameAttribute::NoJoins) &&
            edge_datas.assignment.as_ref().map(|vec| vec.contains(&EdgeAssignment::Join)).unwrap_or(false)
        {
            self.frame_attributes.swap_remove(&FrameAttribute::NoJoins);
        }
        if self.frame_attributes.contains(&FrameAttribute::Manifold) {
            let vh = self.vertices_half_edges.as_ref().unwrap();
            let ec = self.edges_face_corners.as_ref().unwrap();
            for &v in edge_datas.vertices.iter().flatten() {
                if vh[v].iter().all(|&h| ec.at(h).len() + ec.at(h.flipped()).len() == 2) {
                    self.frame_attributes.swap_remove(&FrameAttribute::Orientable);
                    self.frame_attributes.swap_remove(&FrameAttribute::Manifold);
                    break;
                }
            }
        }
        self.add_edges_unchecked(edge_datas)
    }

    /// Gets the data for an edge all in one place.
    pub fn edge_data(&mut self, edge: Edge) -> EdgeData {
        EdgeData {
            vertices: self.edges_vertices.as_ref().unwrap()[edge],
            face_corners: self.edges_face_corners.as_ref().map(|vec| vec[edge].clone()),
            assignment: self.edges_assignment.as_ref().map(|vec| vec[edge].clone()),
            fold_angle_f64: self.edges_fold_angle_f64.as_ref().map(|vec| vec[edge].clone()),
            fold_angle_exact: self.edges_fold_angle_exact.as_ref().map(|vec| vec[edge].clone()),
            length_f64: self.edges_length_f64.as_ref().map(|vec| vec[edge].clone()),
            length2_exact: self.edges_length2_exact.as_ref().map(|vec| vec[edge].clone()),
            custom: self.edges_custom.iter_mut().map(|(k, vec)| (k.clone(), vec[edge].clone())).collect()
        }
    }

    /// Adds an edge with the same data as the edge given, except between `vertices`.
    /// 
    /// Note that `edge_data.face_corners` gets ignored.
    pub fn add_edge_like(&mut self, edge: Edge, vertices: [Vertex; 2]) -> Edge {
        let data = EdgeData { vertices, ..self.edge_data(edge) };
        self.add_edge(data)
    }

    /// Adds a face given some face data and returns its index
    /// without checking whether it preserves the frame attributes.
    /// See [`add_face`](Frame::add_face) for further documentation.
    pub fn add_face_unchecked(&mut self, mut face_data: FaceData) -> Face {
        let index = Face(self.faces_half_edges.as_ref().unwrap().len());

        let ec = self.edges_face_corners.as_mut().unwrap();
        face_data.half_edges.iter().enumerate().for_each(|(i, &h)| ec.at_mut(h).push(FaceCorner(index, i)));
        self.faces_half_edges.as_mut().unwrap().push(face_data.half_edges);
        self.faces_custom.iter_mut().for_each(|(k, vec)| vec.push(face_data.custom.swap_remove(k).unwrap()));
        index
    }

    /// Adds faces given some face data and returns the index of the first face
    /// without checking whether adding the faces preserves the frame attributes.
    /// See [`add_faces`](Frame::add_faces) for further documentation.
    pub fn add_faces_unchecked(&mut self, mut face_datas: FaceDatas) -> Face {
        let index = Face(self.faces_half_edges.as_ref().unwrap().len());
        let _len = face_datas.len();

        let ec = self.edges_face_corners.as_mut().unwrap();
        for (f, hs) in face_datas.half_edges.iter().enumerate() {
            hs.iter().enumerate().for_each(|(i, &h)| ec.at_mut(h).push(FaceCorner(index + Face(f), i)));
        }
        self.faces_half_edges.as_mut().unwrap().extend(face_datas.half_edges);
        self.faces_custom.iter_mut().for_each(|(k, vec)| vec.extend(face_datas.custom.swap_remove(k).unwrap()));
        index
    }

    /// Flips a face. Its half-edges get reversed. 
    pub fn flip_face(&mut self, face: Face) {
        let ec = self.edges_face_corners.as_mut().unwrap();
        let fh = self.faces_half_edges.as_mut().unwrap();

        let num_edges = fh[face].len();
        for (i, h) in fh[face].iter_mut().enumerate() {
            let pos = ec.at(*h).iter().position(|c| *c == FaceCorner(face, i)).unwrap();
            ec.at_mut(*h).remove(pos);
            h.flip();
            ec.at_mut(*h).push(FaceCorner(face, num_edges - i - 1));
        }
        fh[face].reverse();
    }

    fn fix_manifold_attributes_after_adding_faces(&mut self, faces: impl IntoIterator<Item = Face>) {
        let mut fan_map = indexmap! {};

        let fh = self.faces_half_edges.as_ref().unwrap();
        let ev = self.edges_vertices.as_ref().unwrap();
        let mut half_edges = faces.into_iter().flat_map(|f| fh[f].iter().copied()).collect::<Vec<_>>();
        half_edges.sort();
        half_edges.dedup();
        let mut vertices = half_edges.iter().map(|&h| ev.at(h)[0]).collect::<Vec<_>>();
        vertices.sort();
        vertices.dedup();

        if self.frame_attributes.contains(&FrameAttribute::Manifold) {
            for &v in &vertices {
                // For manifold checking,
                // It suffices to check that the vertex is manifold, because if the half-edge isn't manifold,
                // then neither is the vertex.
                // Meanwhile, we also get the fans, since we need them for reordering the half-edges around a vertex.
                let entry = fan_map.entry(v);
                if let Entry::Vacant(_) = &entry {
                    let fans = self.vertex_fans_if_manifold(v);
                    if let Some(fans) = fans {
                        entry.or_insert(fans);
                    } else {
                        self.frame_attributes.swap_remove(&FrameAttribute::Orientable);
                        self.frame_attributes.swap_remove(&FrameAttribute::Manifold);
                        break;
                    }
                }
            }
        }

        if self.frame_attributes.contains(&FrameAttribute::Orientable) {
            let ec = self.edges_face_corners.as_ref().unwrap();
            for &h in &half_edges {
                if ec.at(h).len() == 2 {
                    self.frame_attributes.swap_remove(&FrameAttribute::Orientable);
                    break;
                }
            }
        }

        if self.frame_attributes.contains(&FrameAttribute::Orientable) {
            // Flip fans if necessary
            let ec = self.edges_face_corners.as_ref().unwrap();
            let fh = self.faces_half_edges.as_ref().unwrap();
            for (_, fans) in &mut fan_map {
                for (fan, _) in fans {
                    if fan.len() >= 2 && (ec.at(fan[0]).is_empty() || fh.at(ec.at(fan[0])[0].prev(fh)).flipped() != fan[1]) {
                        fan.reverse();
                    }
                }
            }
        }

        // Finally, resort vertices_half_edges
        if self.frame_attributes.contains(&FrameAttribute::Manifold) {
            let vh = self.vertices_half_edges.as_mut().unwrap();
            for (v, fans) in fan_map {
                vh[v] = fans.into_iter().flat_map(|(fan, _)| fan).collect();
            }
        }
    }

    /// Adds an face given some face data and returns its index.
    /// 
    /// The half-edges in `face_data` should form a loop.
    /// 
    /// If the mesh is no longer a manifold, then [`Manifold`](FrameAttribute::Manifold) is cleared.
    /// If the mesh is no longer orient*ed* (even if it's stil orient*able*), then
    /// [`Orientable`](FrameAttribute::Orientable) gets cleared. You'll have to call
    /// [`try_into_orientable`](Frame::try_into_orientable) again.
    /// 
    /// All fields that exists in this frame must have their respective field specified
    /// in `face_data`
    pub fn add_face(&mut self, face_data: FaceData) -> Face {
        let ev = self.edges_vertices.as_ref().unwrap();
        for (&h0, &h1) in face_data.half_edges.iter().zip(face_data.half_edges.iter().cycle().skip(1)) {
            assert_eq!(ev.at(h0)[1], ev.at(h1)[0], "half-edges in added face must form a loop")
        }
        let index = self.add_face_unchecked(face_data);
        self.fix_manifold_attributes_after_adding_faces(iter::once(index));
        index
    }

    /// Adds faces given some face data and returns the index of the first face.
    /// 
    /// The half-edges in `face_data` should form a loop.
    /// 
    /// If the mesh is no longer a manifold, then [`Manifold`](FrameAttribute::Manifold) is cleared.
    /// If the mesh is no longer orient*ed* (even if it's stil orient*able*), then
    /// [`Orientable`](FrameAttribute::Orientable) gets cleared. You'll have to call
    /// [`try_into_orientable`](Frame::try_into_orientable) again.
    /// 
    /// All fields that exists in this frame must have their respective field specified
    /// in `face_data`
    pub fn add_faces(&mut self, face_datas: FaceDatas) -> Face {
        let ev = self.edges_vertices.as_ref().unwrap();
        for hs in &face_datas.half_edges {
            for (&h0, &h1) in hs.iter().zip(hs.iter().cycle().skip(1)) {
                assert_eq!(ev.at(h0)[1], ev.at(h1)[0], "half-edges in added face must form a loop")
            }
        }
        let len = face_datas.len();
        let index = self.add_faces_unchecked(face_datas);
        self.fix_manifold_attributes_after_adding_faces((0..len).map(|i| index + Face(i)));
        index
    }

    ///// Tries to add a face given some face data and vertices (instead of half-edges) and returns its index.
    ///// Modifies frame attributes as necessary to remove contradictions.
    ///// Returns `None` if the resulting face is ambiguous (multiple edges connect the same two adjacent vertices.)
    ///// 
    ///// The half-edges in `face_data` should form a loop.
    ///// 
    ///// All fields that exists in this frame must have their respective field specified
    ///// in `face_data`
    //pub fn add_face_from_vertices(&mut self, vertices: &[Vertex], other_data: FaceData) -> Option<Face> {
    //    
    //}

    /// Gets the data for a face all in one place.
    pub fn face_data(&mut self, face: Face) -> FaceData {
        FaceData {
            half_edges: self.faces_half_edges.as_ref().unwrap()[face].clone(),
            custom: self.faces_custom.iter_mut().map(|(k, vec)| (k.clone(), vec[face].clone())).collect()
        }
    }

    /// Adds a face with the same data as the face given, except from `half_edges`.
    /// 
    /// Note that `face_data.face_corners` gets ignored.
    pub fn add_face_like(&mut self, face: Face, half_edges: Vec<HalfEdge>) -> Face {
        let data = FaceData { half_edges, ..self.face_data(face) };
        self.add_face(data)
    }

    fn remap_vertices_refs(&mut self, new_indices: TiVec<Vertex, Option<Vertex>>) {
        if let Some(edges_vertices) = self.edges_vertices.as_mut() {
            for (e, vertices) in edges_vertices.iter_mut_enumerated() {
                *vertices = vertices
                    .map(|v| new_indices[v].unwrap_or_else(|| panic!("tried to remap vertex {v} to None while an edge {e} contains it")))
            }
        }
    }

    fn remap_edges_refs(&mut self, new_indices: TiVec<Edge, Option<Edge>>) {
        if let Some(faces_half_edges) = self.faces_half_edges.as_mut() {
            for (f, half_edges) in faces_half_edges.iter_mut_enumerated() {
                *half_edges = half_edges.drain(..)
                    .map(|h| new_indices[h.edge()].map(|e| HalfEdge::new(e, h.flip_bit()))
                        .unwrap_or_else(|| panic!("tried to remap half-edge {h} to None while a face {f} contains it")))
                    .collect::<Vec<_>>();
            }
        }

        if let Some(vertices_half_edges) = self.vertices_half_edges.as_mut() {
            for half_edges in vertices_half_edges {
                *half_edges = half_edges.drain(..)
                    .flat_map(|h| new_indices[h.edge()].map(|e| HalfEdge::new(e, h.flip_bit())))
                    .collect::<Vec<_>>();
            }
        }

        if let Some(edge_orders) = self.edge_orders.as_mut() {
            *edge_orders = edge_orders.drain(..)
                .flat_map(|(e1, e2, order)|
                    new_indices[e1].and_then(|e1| new_indices[e2].map(|e2| (e1, e2, order))))
                .collect::<Vec<_>>();
        }
    }

    fn remap_faces_refs(&mut self, new_indices: TiVec<Face, Option<Face>>) {
        if let Some(edges_face_corners) = self.edges_face_corners.as_mut() {
            for (_, corners) in edges_face_corners.half_iter_mut_enumerated() {
                *corners = corners.drain(..)
                    .flat_map(|c| new_indices[c.face()].map(|f| FaceCorner(f, c.corner())))
                    .collect::<Vec<_>>();
            }
        }

        if let Some(face_orders) = self.face_orders.as_mut() {
            *face_orders = face_orders.drain(..)
                .flat_map(|(f1, f2, order)|
                    new_indices[f1].and_then(|f1| new_indices[f2].map(|f2| (f1, f2, order))))
                .collect::<Vec<_>>();
        }
    }

    /// Removes all edges and faces attached to the vertex as well and returns their data.
    fn remove_vertex_generic(
        &mut self,
        vertex: Vertex,
        strategy: impl RemoveStrategy,
    ) -> (VertexData, IndexMap<Edge, EdgeData>, IndexMap<Face, FaceData>) {
        let face_data = if let Some(edges_face_corners) = self.edges_face_corners.as_ref() {
            let mut faces = self.vertices_half_edges.as_ref().unwrap()[vertex].iter()
                .flat_map(|&h| edges_face_corners.at(h).iter().map(|c| c.face()))
                .collect::<Vec<_>>();
            faces.sort();
            faces.dedup();
            faces.reverse();
            faces.into_iter().map(|f| (f, self.remove_face_generic(f, strategy.clone()))).collect::<IndexMap<_, _>>()
        } else { indexmap! {} };

        let edge_data = if let Some(vertices_half_edges) = self.vertices_half_edges.as_ref() {
            let mut edges = vertices_half_edges[vertex].iter().map(|h| h.edge()).collect::<Vec<_>>();
            edges.sort();
            edges.reverse(); // no need for dedup because self-loops aren't allowed
            edges.into_iter().map(|e| (e, self.remove_edge_generic(e, strategy.clone()).0)).collect::<IndexMap<_, _>>()
        } else { indexmap! {} };

        let new_indices = (0..self.num_vertices).map(
            |v| strategy.index_map(self.num_vertices, vertex, Vertex(v))).collect::<TiVec<_, _>>();
        self.num_vertices -= 1;
        
        let data = VertexData {
            coords_f64: self.vertices_coords_f64.as_mut().map(|vec| strategy.remove_column(vec, vertex.0)),
            coords_exact: self.vertices_coords_exact.as_mut().map(|vec| strategy.remove_column(vec, vertex.0)),
            half_edges: self.vertices_half_edges.as_mut().map(|vec| strategy.remove(vec, vertex)),
            custom: self.vertices_custom.iter_mut()
                .map(|(field, value)| (field.clone(), strategy.remove(value, vertex)))
                    .collect::<IndexMap<_, _>>()
        };
        self.remap_vertices_refs(new_indices);
        
        (data, edge_data, face_data)
    }

    /// Shift-removes a vertex `vertex`, shifting all vertexs after it to fill the gap,
    /// and shift-removes all edges and faces attached to the vertex.
    /// Returns the vertex data for that vertex, and the edge and face data for all removed edges and faces.
    pub fn shift_remove_vertex(&mut self, vertex: Vertex) -> (VertexData, IndexMap<Edge, EdgeData>, IndexMap<Face, FaceData>) {
        self.remove_vertex_generic(vertex, ShiftRemove)
    }

    /// Swap-removes a vertex `vertex`, moving the last vertex into its position,
    /// and swap-removes all edges and faces attached to the vertex.
    /// Returns the vertex data for that vertex, and the edge and face data for all removed edges and faces.
    pub fn swap_remove_vertex(&mut self, vertex: Vertex) -> (VertexData, IndexMap<Edge, EdgeData>, IndexMap<Face, FaceData>) {
        self.remove_vertex_generic(vertex, SwapRemove)
    }

    /// Removes all faces attached to the edge as well and returns their data.
    fn remove_edge_generic(
        &mut self,
        edge: Edge,
        strategy: impl RemoveStrategy,
    ) -> (EdgeData, IndexMap<Face, FaceData>) {
        let face_data = if let Some(edges_face_corners) = self.edges_face_corners.as_ref() {
            let mut faces = edges_face_corners[edge].iter().flatten().map(|c| c.face()).collect::<Vec<_>>();
            faces.sort();
            faces.dedup();
            faces.reverse();
            faces.into_iter().map(|f| (f, self.remove_face_generic(f, strategy.clone()))).collect::<IndexMap<_, _>>()
        } else { indexmap! {} };

        let edges_vertices = self.edges_vertices.as_mut().expect("tried to remove an edge without there being any");

        let new_indices = edges_vertices.iter_enumerated().map(
            |(e, _)| strategy.index_map(edges_vertices.len(), edge, e)).collect::<TiVec<_, _>>();
        
        let data = EdgeData {
            vertices: strategy.remove(edges_vertices, edge),
            face_corners: self.edges_face_corners.as_mut().map(|vec| strategy.remove(vec, edge)),
            assignment: self.edges_assignment.as_mut().map(|vec| strategy.remove(vec, edge)),
            fold_angle_f64: self.edges_fold_angle_f64.as_mut().map(|vec| strategy.remove(vec, edge)),
            fold_angle_exact: self.edges_fold_angle_exact.as_mut().map(|vec| strategy.remove(vec, edge)),
            length_f64: self.edges_length_f64.as_mut().map(|vec| strategy.remove(vec, edge)),
            length2_exact: self.edges_length2_exact.as_mut().map(|vec| strategy.remove(vec, edge)),
            custom: self.edges_custom.iter_mut()
                .map(|(field, value)| (field.clone(), strategy.remove(value, edge)))
                    .collect::<IndexMap<_, _>>()
        };
        self.remap_edges_refs(new_indices);
        
        (data, face_data)
    }

    /// Shift-removes a edge `edge`, shifting all edges after it to fill the gap,
    /// and shift-removes all faces attached to the edge.
    /// Returns the edge data for that edge, and the face data for all removed faces.
    pub fn shift_remove_edge(&mut self, edge: Edge) -> (EdgeData, IndexMap<Face, FaceData>) {
        self.remove_edge_generic(edge, ShiftRemove)
    }

    /// Swap-removes a edge `edge`, moving the last edge into its position,
    /// and swap-removes all faces attached to the edge.
    /// Returns the edge data for that edge, and the face data for all removed faces.
    pub fn swap_remove_edge(&mut self, edge: Edge) -> (EdgeData, IndexMap<Face, FaceData>) {
        self.remove_edge_generic(edge, SwapRemove)
    }

    fn remove_face_generic(
        &mut self,
        face: Face,
        strategy: impl RemoveStrategy,
    ) -> FaceData {
        let faces_half_edges = self.faces_half_edges.as_mut().expect("tried to remove a face without there being any");

        let new_indices = faces_half_edges.iter_enumerated().map(
            |(f, _)| strategy.index_map(faces_half_edges.len(), face, f)).collect::<TiVec<_, _>>();
        
        let half_edges = strategy.remove(faces_half_edges, face);
        let custom = self.faces_custom.iter_mut()
            .map(|(field, value)| (field.clone(), strategy.remove(value, face)))
            .collect::<IndexMap<_, _>>();
        self.remap_faces_refs(new_indices);
        
        FaceData { half_edges, custom }
    }

    /// Shift-removes a face `face`, shifting all faces after it to fill the gap.
    /// Returns the face data for that face.
    pub fn shift_remove_face(&mut self, face: Face) -> FaceData {
        self.remove_face_generic(face, ShiftRemove)
    }

    /// Swap-removes a face `face`, moving the last face into its position.
    /// Returns the face data for that face.
    pub fn swap_remove_face(&mut self, face: Face) -> FaceData {
        self.remove_face_generic(face, SwapRemove)
    }
}

/// Given two edges defined by their vertices, gets a vertex incident to both of them, if one exists.
pub fn edges_vertices_incident(e1: [Vertex; 2], e2: [Vertex; 2]) -> Option<Vertex> {
    e1.into_iter().flat_map(|v1| e2.into_iter().find(|v2| v1 == *v2)).next()
}

/// Given two vertices defined by their half-edges, gets the half-edges incident to both of them.
/// The half-edges point in the v1->v2 direction.
/// 
/// Make sure the vertices aren't the exact same vertex before calling this.
pub fn vertices_half_edges_incident(v1: &[HalfEdge], v2: &[HalfEdge]) -> Vec<HalfEdge> {
    v1.iter().copied().filter(|h1| v2.iter().find(|h2| h2.edge() == h1.edge()).is_some()).collect::<Vec<_>>()
}

/// Given an edge and a vertex on that edge, gets the other vertex
pub fn try_other_vertex(edge: [Vertex; 2], vertex: Vertex) -> Option<Vertex> {
    if edge[0] == vertex { Some(edge[1]) } else if edge[1] == vertex { Some(edge[0]) } else { None }
}

/// Given an edge and a vertex on that edge, gets the other vertex
pub fn other_vertex(edge: [Vertex; 2], vertex: Vertex) -> Vertex {
    if edge[0] == vertex { edge[1] } else { edge[0] }
}

#[cfg(test)]
mod test {
    use exact_number::{based_expr, Angle};
    use indexmap::{indexmap, indexset};
    use nalgebra::{DMatrix, DVector};
    use serde_json::json;
    use typed_index_collections::ti_vec;

    use crate::{filter::{edges_vertices_incident, other_vertex}, fold::{EdgeAssignment, EdgeData, EdgeOrder, FaceData, FaceOrder, Frame, FrameAttribute, VertexData}};
    use crate::fold::{Vertex as V, Edge as E, Face as F, HalfEdge as H, FaceCorner as C};

    #[test]
    fn test_edges_vertices_incident() {
        assert_eq!(edges_vertices_incident([V(0), V(2)], [V(2), V(4)]), Some(V(2)));
        assert_eq!(edges_vertices_incident([V(0), V(2)], [V(4), V(2)]), Some(V(2)));
        assert_eq!(edges_vertices_incident([V(2), V(0)], [V(2), V(4)]), Some(V(2)));
        assert_eq!(edges_vertices_incident([V(2), V(0)], [V(4), V(2)]), Some(V(2)));
        assert_eq!(edges_vertices_incident([V(2), V(0)], [V(4), V(5)]), None);
    }

    #[test]
    fn test_other_vertex() {
        assert_eq!(other_vertex([V(5), V(4)], V(4)), V(5));
        assert_eq!(other_vertex([V(5), V(4)], V(5)), V(4));
        assert_eq!(other_vertex([V(3), V(7)], V(3)), V(7));
        assert_eq!(other_vertex([V(3), V(7)], V(7)), V(3));
    }

    #[test]
    fn test_add_vertex() {
        let mut frame = Frame {
            num_vertices: 1,
            vertices_coords_f64: Some(DMatrix::from_vec(2, 1, vec![
                0.0, 1.0,
            ])),
            vertices_coords_exact: Some(DMatrix::from_vec(2, 1, vec![
                based_expr!(0), based_expr!(1),
            ])),
            vertices_half_edges: Some(ti_vec![vec![]]),
            vertices_custom: indexmap! {
                "test:test".to_owned() => ti_vec![json!(1)]
            },
            ..Default::default()
        };
        let index = frame.add_vertex(VertexData {
            coords_f64: Some(DVector::from_vec(vec![2.0, 1.0])),
            coords_exact: Some(DVector::from_vec(vec![based_expr!(2), based_expr!(1)])),
            half_edges: None,
            custom: indexmap! { "test:test".to_owned() => json!(3) },
        });
        assert_eq!(index, V(1));
        assert_eq!(frame.num_vertices, 2);
        assert_eq!(frame.vertices_coords_f64, Some(DMatrix::from_vec(2, 2, vec![
            0.0, 1.0,
            2.0, 1.0,
        ])));
        assert_eq!(frame.vertices_coords_exact, Some(DMatrix::from_vec(2, 2, vec![
            based_expr!(0), based_expr!(1),
            based_expr!(2), based_expr!(1),
        ])));
        assert_eq!(frame.vertices_half_edges, Some(ti_vec![vec![], vec![]]));
        assert_eq!(frame.vertices_custom, indexmap! {
            "test:test".to_owned() => ti_vec![json!(1), json!(3)]
        });
    }

    #[test]
    fn test_add_edge() {
        use EdgeAssignment::*;

        let mut frame = Frame {
            num_vertices: 3,
            edges_assignment: Some(ti_vec![Boundary, Boundary]),
            edges_fold_angle_f64: Some(ti_vec![0.0, 0.0]),
            edges_fold_angle_exact: Some(ti_vec![Angle::ZERO, Angle::ZERO]),
            edges_length_f64: Some(ti_vec![2.0, 2.5]),
            edges_length2_exact: Some(ti_vec![based_expr!(4), based_expr!(25/4)]),
            edges_custom: indexmap! {
                "test:test".to_owned() => ti_vec![json!(1), json!(2)]
            },
            ..Default::default()
        }.with_topology_vh_fh(Some(ti_vec![
            vec![H(0)],
            vec![H(1), H(2)],
            vec![H(3)],
        ]), None);
        let index = frame.add_edge(EdgeData {
            vertices: [V(2), V(0)],
            face_corners: None,
            assignment: Some(EdgeAssignment::Unassigned),
            fold_angle_f64: Some(0.0),
            fold_angle_exact: Some(Angle::ZERO),
            length_f64: Some(1.0),
            length2_exact: Some(based_expr!(1)),
            custom: indexmap! { "test:test".to_owned() => json!(3) },
        });
        assert_eq!(index, E(2));
        frame.vertices_half_edges.as_mut().unwrap().iter_mut().for_each(|hs| hs.sort());
        assert_eq!(frame.vertices_half_edges, Some(ti_vec![
            vec![H(0), H(5)],
            vec![H(1), H(2)],
            vec![H(3), H(4)],
        ]));
        assert_eq!(frame.edges_vertices, Some(ti_vec![
            [V(0), V(1)],
            [V(1), V(2)],
            [V(2), V(0)],
        ]));
        assert_eq!(frame.edges_face_corners, None);
        assert_eq!(frame.edges_assignment, Some(ti_vec![Boundary, Boundary, Unassigned]));
        assert_eq!(frame.edges_fold_angle_f64, Some(ti_vec![0.0, 0.0, 0.0]));
        assert_eq!(frame.edges_fold_angle_exact, Some(ti_vec![Angle::ZERO, Angle::ZERO, Angle::ZERO]));
        assert_eq!(frame.edges_length_f64, Some(ti_vec![2.0, 2.5, 1.0]));
        assert_eq!(frame.edges_length2_exact, Some(ti_vec![based_expr!(4), based_expr!(25/4), based_expr!(1)]));
        assert_eq!(frame.edges_custom, indexmap! {
            "test:test".to_owned() => ti_vec![json!(1), json!(2), json!(3)]
        });

        let mut frame = Frame {
            num_vertices: 3,
            edges_assignment: Some(ti_vec![Boundary, Boundary]),
            frame_attributes: indexset! {FrameAttribute::NoCuts, FrameAttribute::NoJoins},
            ..Default::default()
        }.with_topology_vh_fh(Some(ti_vec![
            vec![H(0)],
            vec![H(1), H(2)],
            vec![H(3)],
        ]), None);
        frame.add_edge(EdgeData { assignment: Some(Cut), ..EdgeData::default_with_vertices([V(2), V(0)]) });
        assert_eq!(frame.frame_attributes, indexset! {FrameAttribute::NoJoins});
        frame.add_edge(EdgeData { assignment: Some(Join), ..EdgeData::default_with_vertices([V(2), V(0)]) });
        assert_eq!(frame.frame_attributes, indexset! {});

        let mut frame = Frame {
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
        ])).try_into_orientable().unwrap().0;
        frame.add_edge(EdgeData { assignment: Some(Cut), ..EdgeData::default_with_vertices([V(0), V(2)]) });
        assert_eq!(frame.frame_attributes, indexset! {FrameAttribute::Manifold, FrameAttribute::Orientable});
        frame.add_edge(EdgeData { assignment: Some(Join), ..EdgeData::default_with_vertices([V(4), V(0)]) });
        assert_eq!(frame.frame_attributes, indexset! {});
    }

    #[test]
    fn test_add_face() {
        use FrameAttribute::*;

        // 0---3---1---4---2
        // |       |       |
        // 0       1   0   2
        // |       |       |
        // 3---5---4---6---5
        // all edges point down and right
        // all faces start at the top-left corner and go ccw
        let mut frame = Frame {
            faces_custom: indexmap! {
                "test:test".to_owned() => ti_vec![json!(1)]
            },
            ..Default::default()
        }.with_topology_vh_fh(Some(ti_vec![
            vec![H(0), H(6)],
            vec![H(7), H(2), H(8)],
            vec![H(9), H(4)],
            vec![H(10), H(1)],
            vec![H(12), H(3), H(11)],
            vec![H(5), H(13)],
        ]), Some(ti_vec![
            vec![H(2), H(12), H(5), H(9)],
        ])).try_into_orientable().unwrap().0;
        
        let index = frame.add_face(FaceData {
            half_edges: vec![H(0), H(10), H(3), H(7)],
            custom: indexmap! { "test:test".to_owned() => json!(3) }
        });
        assert_eq!(frame.frame_attributes, indexset! {Manifold, Orientable});
        frame.assert_topology_consistent();
        assert_eq!(index, F(1));
        assert_eq!(frame.faces_half_edges, Some(ti_vec![
            vec![H(2), H(12), H(5), H(9)],
            vec![H(0), H(10), H(3), H(7)]
        ]));
        assert_eq!(frame.faces_custom, indexmap! {
            "test:test".to_owned() => ti_vec![json!(1), json!(3)]
        });
        let two_panel = frame;

        let mut frame = two_panel.clone();
        // Still orientable, but since it needs a flip, the flag gets cleared.
        let [h14, _h15] = frame.add_edge_like(E(0), [V(2), V(0)]).split();
        let [_h16, h17] = frame.add_edge_like(E(0), [V(5), V(3)]).split();
        let index = frame.add_face(FaceData {
            half_edges: vec![H(5), h14, H(0), h17],
            custom: indexmap! { "test:test".to_owned() => json!(3) }
        });
        assert_eq!(frame.frame_attributes, indexset! {Manifold});
        frame.assert_topology_consistent();
        assert_eq!(index, F(2));
        assert_eq!(frame.faces_half_edges, Some(ti_vec![
            vec![H(2), H(12), H(5), H(9)],
            vec![H(0), H(10), H(3), H(7)],
            vec![H(5), h14, H(0), h17],
        ]));

        let mut frame = two_panel.clone();
        // Not orientable, but is still a manifold
        let [h14, _h15] = frame.add_edge_like(E(0), [V(2), V(3)]).split();
        let [_h16, h17] = frame.add_edge_like(E(0), [V(5), V(0)]).split();
        let index = frame.add_face(FaceData {
            half_edges: vec![H(5), h14, H(1), h17],
            custom: indexmap! { "test:test".to_owned() => json!(3) }
        });
        assert_eq!(frame.frame_attributes, indexset! {Manifold});
        frame.assert_topology_consistent();
        assert_eq!(index, F(2));
        assert_eq!(frame.faces_half_edges, Some(ti_vec![
            vec![H(2), H(12), H(5), H(9)],
            vec![H(0), H(10), H(3), H(7)],
            vec![H(5), h14, H(1), h17],
        ]));

        let mut frame = two_panel.clone();
        // Not even a manifold
        let index = frame.add_face(FaceData {
            half_edges: vec![H(0), H(10), H(3), H(7)],
            custom: indexmap! { "test:test".to_owned() => json!(3) }
        });
        assert_eq!(frame.frame_attributes, indexset! {});
        frame.assert_topology_consistent();
        assert_eq!(index, F(2));
        assert_eq!(frame.faces_half_edges, Some(ti_vec![
            vec![H(2), H(12), H(5), H(9)],
            vec![H(0), H(10), H(3), H(7)],
            vec![H(0), H(10), H(3), H(7)],
        ]));
    }

    #[test]
    fn test_swap_remove_vertex() {
        // TODO: Add other test cases, and test shift_remove_vertex. I'm getting tired
        // of how long it takes to write these test cases...

        // 0---3---1---4---2
        // |       |       |
        // 0   0   1   1   2
        // |       |       |
        // 3---5---4---6---5
        // all edges point down and right
        // all faces start at the top-left corner and go ccw
        let wrap_fold = Frame {
            num_vertices: 6,
            vertices_coords_f64: Some(DMatrix::from_vec(2, 6, vec![
                0.0, 1.0,
                2.0, 1.0,
                4.0, 1.0,
                0.0, 0.0,
                2.0, 0.0,
                4.0, 0.0,
            ])),
            vertices_coords_exact: Some(DMatrix::from_vec(2, 6, vec![
                based_expr!(0), based_expr!(1),
                based_expr!(2), based_expr!(1),
                based_expr!(4), based_expr!(1),
                based_expr!(0), based_expr!(0),
                based_expr!(2), based_expr!(0),
                based_expr!(4), based_expr!(0),
            ])),
            vertices_custom: indexmap! {
                "test:test".to_owned() => ti_vec![json!(1), json!(2), json!(3), json!(4), json!(5), json!(6)]
            },
            ..Default::default()
        }.with_topology_vh_fh(Some(ti_vec![
            vec![H(0), H(6)],
            vec![H(7), H(2), H(8)],
            vec![H(9), H(4)],
            vec![H(10), H(1)],
            vec![H(12), H(3), H(11)],
            vec![H(5), H(13)],
        ]), Some(ti_vec![
            vec![H(0), H(10), H(3), H(7)],
            vec![H(2), H(12), H(5), H(9)],
        ]));

        let mut frame = wrap_fold.clone();
        let (vertex_data, edge_datas, face_datas) = frame.swap_remove_vertex(V(0));
        assert_eq!(vertex_data, VertexData {
            coords_f64: Some(DVector::from_vec(vec![0.0, 1.0])),
            coords_exact: Some(DVector::from_vec(vec![based_expr!(0), based_expr!(1)])),
            half_edges: Some(vec![]),
            custom: indexmap! { "test:test".to_owned() => json!(1) }
        });
        assert_eq!(edge_datas, indexmap! {
            E(0) => EdgeData {
                vertices: [V(0), V(3)],
                face_corners: Some([vec![], vec![]]),
                assignment: None,
                fold_angle_f64: None,
                fold_angle_exact: None,
                length_f64: None,
                length2_exact: None,
                custom: indexmap! {}
            },
            E(3) => EdgeData {
                vertices: [V(0), V(1)],
                face_corners: Some([vec![], vec![]]),
                assignment: None,
                fold_angle_f64: None,
                fold_angle_exact: None,
                length_f64: None,
                length2_exact: None,
                custom: indexmap! {}
            },
        });
        assert_eq!(face_datas, indexmap! {
            F(0) => FaceData {
                half_edges: vec![H(0), H(10), H(3), H(7)],
                custom: indexmap! {}
            }
        });
        assert_eq!(frame.num_vertices, 5);
        assert_eq!(frame.edges_vertices, Some(ti_vec![
            [V(3), V(4)],
            [V(1), V(4)],
            [V(2), V(0)],
            [V(4), V(0)],
            [V(1), V(2)],
        ]));
        assert_eq!(frame.vertices_half_edges, Some(ti_vec![
            vec![H(5), H(7)],
            vec![H(2), H(8)],
            vec![H(9), H(4)],
            vec![H(0)],
            vec![H(6), H(3), H(1)],
        ]));
        assert_eq!(frame.vertices_coords_f64, Some(DMatrix::from_vec(2, 5, vec![
            4.0, 0.0,
            2.0, 1.0,
            4.0, 1.0,
            0.0, 0.0,
            2.0, 0.0,
        ])));
        assert_eq!(frame.vertices_coords_exact, Some(DMatrix::from_vec(2, 5, vec![
            based_expr!(4), based_expr!(0),
            based_expr!(2), based_expr!(1),
            based_expr!(4), based_expr!(1),
            based_expr!(0), based_expr!(0),
            based_expr!(2), based_expr!(0),
        ])));
        assert_eq!(frame.vertices_custom, indexmap! {
            "test:test".to_owned() => ti_vec![json!(6), json!(2), json!(3), json!(4), json!(5)]
        });
    }

    #[test]
    fn test_shift_remove_edge() {
        use EdgeAssignment::*;
        
        // 0---3---1---4---2
        // |       |       |
        // 0   0   1   1   2
        // |       |       |
        // 3---5---4---6---5
        // all edges point down and right
        // all faces start at the top-left corner and go ccw
        let wrap_fold = Frame {
            edges_assignment: Some(ti_vec![Boundary, Valley, Boundary, Boundary, Boundary, Boundary, Boundary]),
            edges_fold_angle_f64: Some(ti_vec![0.0, -180.0, 0.0, 0.0, 0.0, 0.0, 0.0]),
            edges_fold_angle_exact: Some(ti_vec![
                Angle::ZERO, -Angle::DEG_180, Angle::ZERO, Angle::ZERO, Angle::ZERO, Angle::ZERO, Angle::ZERO,
            ]),
            edges_length_f64: Some(ti_vec![1.0, 1.0, 1.0, 2.0, 2.0, 2.0, 2.0]),
            edges_length2_exact: Some(ti_vec![
                based_expr!(1), based_expr!(1), based_expr!(1), based_expr!(2), based_expr!(2), based_expr!(2), based_expr!(2)
            ]),
            edge_orders: Some(vec![
                (E(0), E(1), EdgeOrder::Right),
                (E(1), E(0), EdgeOrder::Left),
                (E(1), E(2), EdgeOrder::Right),
                (E(2), E(1), EdgeOrder::Left),
            ]),
            edges_custom: indexmap! {
                "test:test".to_owned() => ti_vec![json!(1), json!(2), json!(3), json!(4), json!(5), json!(6), json!(7)]
            },
            ..Default::default()
        }.with_topology_vh_fh(Some(ti_vec![
            vec![H(0), H(6)],
            vec![H(7), H(2), H(8)],
            vec![H(9), H(4)],
            vec![H(10), H(1)],
            vec![H(12), H(3), H(11)],
            vec![H(5), H(13)],
        ]), Some(ti_vec![
            vec![H(0), H(10), H(3), H(7)],
            vec![H(2), H(12), H(5), H(9)],
        ]));

        let mut frame = wrap_fold.clone();
        let (edge_data, face_datas) = frame.shift_remove_edge(E(0));
        assert_eq!(edge_data, EdgeData {
            vertices: [V(0), V(3)],
            face_corners: Some([vec![], vec![]]),
            assignment: Some(EdgeAssignment::Boundary),
            fold_angle_f64: Some(0.0),
            fold_angle_exact: Some(Angle::ZERO),
            length_f64: Some(1.0),
            length2_exact: Some(based_expr!(1)),
            custom: indexmap! { "test:test".to_owned() => json!(1) }
        });
        assert_eq!(face_datas, indexmap! {
            F(0) => FaceData {
                half_edges: vec![H(0), H(10), H(3), H(7)],
                custom: indexmap! {}
            }
        });
        assert_eq!(frame.vertices_half_edges, Some(ti_vec![
            vec![H(4)],
            vec![H(5), H(0), H(6)],
            vec![H(7), H(2)],
            vec![H(8)],
            vec![H(10), H(1), H(9)],
            vec![H(3), H(11)],
        ]));
        assert_eq!(frame.faces_half_edges, Some(ti_vec![
            vec![H(0), H(10), H(3), H(7)],
        ]));
        assert_eq!(frame.edges_vertices, Some(ti_vec![
            [V(1), V(4)],
            [V(2), V(5)],
            [V(0), V(1)],
            [V(1), V(2)],
            [V(3), V(4)],
            [V(4), V(5)],
        ]));
        assert_eq!(frame.edges_face_corners, Some(ti_vec![
            [vec![C(F(0), 0)], vec![]],
            [vec![], vec![C(F(0), 2)]],
            [vec![], vec![]],
            [vec![], vec![C(F(0), 3)]],
            [vec![], vec![]],
            [vec![C(F(0), 1)], vec![]],
        ]));
        assert_eq!(frame.edges_assignment, Some(ti_vec![Valley, Boundary, Boundary, Boundary, Boundary, Boundary]));
        assert_eq!(frame.edges_fold_angle_f64, Some(ti_vec![-180.0, 0.0, 0.0, 0.0, 0.0, 0.0]));
        assert_eq!(frame.edges_fold_angle_exact, Some(ti_vec![
            -Angle::DEG_180, Angle::ZERO, Angle::ZERO, Angle::ZERO, Angle::ZERO, Angle::ZERO,
        ]));
        assert_eq!(frame.edges_length_f64, Some(ti_vec![1.0, 1.0, 2.0, 2.0, 2.0, 2.0]));
        assert_eq!(frame.edges_length2_exact, Some(ti_vec![
            based_expr!(1), based_expr!(1), based_expr!(2), based_expr!(2), based_expr!(2), based_expr!(2)
        ]));
        assert_eq!(frame.edge_orders, Some(vec![
            (E(0), E(1), EdgeOrder::Right),
            (E(1), E(0), EdgeOrder::Left),
        ]));
        assert_eq!(frame.edges_custom, indexmap! {
            "test:test".to_owned() => ti_vec![json!(2), json!(3), json!(4), json!(5), json!(6), json!(7)]
        });

        // Remove the edge adjacent to both faces
        let mut frame = wrap_fold.clone();
        let (edge_data, face_datas) = frame.shift_remove_edge(E(1));
        assert_eq!(edge_data, EdgeData {
            vertices: [V(1), V(4)],
            face_corners: Some([vec![], vec![]]),
            assignment: Some(EdgeAssignment::Valley),
            fold_angle_f64: Some(-180.0),
            fold_angle_exact: Some(-Angle::DEG_180),
            length_f64: Some(1.0),
            length2_exact: Some(based_expr!(1)),
            custom: indexmap! { "test:test".to_owned() => json!(2) }
        });
        assert_eq!(face_datas, indexmap! {
            F(0) => FaceData {
                half_edges: vec![H(0), H(10), H(3), H(7)],
                custom: indexmap! {}
            },
            F(1) => FaceData {
                half_edges: vec![H(2), H(12), H(5), H(9)],
                custom: indexmap! {}
            }
        });
        assert_eq!(frame.vertices_half_edges, Some(ti_vec![
            vec![H(0), H(4)],
            vec![H(5), H(6)],
            vec![H(7), H(2)],
            vec![H(8), H(1)],
            vec![H(10), H(9)],
            vec![H(3), H(11)],
        ]));
        assert_eq!(frame.faces_half_edges, Some(ti_vec![]));
        assert_eq!(frame.edges_vertices, Some(ti_vec![
            [V(0), V(3)],
            [V(2), V(5)],
            [V(0), V(1)],
            [V(1), V(2)],
            [V(3), V(4)],
            [V(4), V(5)],
        ]));
        assert_eq!(frame.edges_face_corners, Some(ti_vec![
            [vec![], vec![]],
            [vec![], vec![]],
            [vec![], vec![]],
            [vec![], vec![]],
            [vec![], vec![]],
            [vec![], vec![]],
        ]));
        assert_eq!(frame.edges_assignment, Some(ti_vec![Boundary, Boundary, Boundary, Boundary, Boundary, Boundary]));
        assert_eq!(frame.edges_fold_angle_f64, Some(ti_vec![0.0, 0.0, 0.0, 0.0, 0.0, 0.0]));
        assert_eq!(frame.edges_fold_angle_exact, Some(ti_vec![
            Angle::ZERO, Angle::ZERO, Angle::ZERO, Angle::ZERO, Angle::ZERO, Angle::ZERO
        ]));
        assert_eq!(frame.edges_length_f64, Some(ti_vec![1.0, 1.0, 2.0, 2.0, 2.0, 2.0]));
        assert_eq!(frame.edges_length2_exact, Some(ti_vec![
            based_expr!(1), based_expr!(1), based_expr!(2), based_expr!(2), based_expr!(2), based_expr!(2)
        ]));
        assert_eq!(frame.edge_orders, Some(vec![]));
        assert_eq!(frame.edges_custom, indexmap! {
            "test:test".to_owned() => ti_vec![json!(1), json!(3), json!(4), json!(5), json!(6), json!(7)]
        });

        // Slitted face, remove edge adjacent to the same face multiple times
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
        let (edge_data, face_datas) = frame.shift_remove_edge(E(8));
        assert_eq!(edge_data, EdgeData {
            vertices: [V(0), V(4)],
            face_corners: Some([vec![], vec![]]),
            assignment: None,
            fold_angle_f64: None,
            fold_angle_exact: None,
            length_f64: None,
            length2_exact: None,
            custom: indexmap! {}
        });
        assert_eq!(face_datas, indexmap! {
            F(0) => FaceData {
                half_edges: vec![H(0), H(2), H(4), H(6), H(16), H(15), H(13), H(11), H(9), H(17)],
                custom: indexmap! {}
            },
        });
        assert_eq!(frame.vertices_half_edges, Some(ti_vec![
            vec![H(0), H(7)],
            vec![H(1), H(2)],
            vec![H(3), H(4)],
            vec![H(5), H(6)],
            vec![H(8), H(15)],
            vec![H(9), H(10)],
            vec![H(11), H(12)],
            vec![H(13), H(14)],
        ]));
        assert_eq!(frame.faces_half_edges, Some(ti_vec![
            vec![H(8), H(10), H(12), H(14)],
        ]));
        assert_eq!(frame.edges_vertices, Some(ti_vec![
            [V(0), V(1)],
            [V(1), V(2)],
            [V(2), V(3)],
            [V(3), V(0)],
            [V(4), V(5)],
            [V(5), V(6)],
            [V(6), V(7)],
            [V(7), V(4)],
        ]));
        assert_eq!(frame.edges_face_corners, Some(ti_vec![
            [vec![], vec![]],
            [vec![], vec![]],
            [vec![], vec![]],
            [vec![], vec![]],
            [vec![C(F(0), 0)], vec![]],
            [vec![C(F(0), 1)], vec![]],
            [vec![C(F(0), 2)], vec![]],
            [vec![C(F(0), 3)], vec![]],
        ]));
    }

    #[test]
    fn test_swap_remove_edge() {
        use EdgeAssignment::*;
        
        // 0---3---1---4---2
        // |       |       |
        // 0   0   1   1   2
        // |       |       |
        // 3---5---4---6---5
        // all edges point down and right
        // all faces start at the top-left corner and go ccw
        let wrap_fold = Frame {
            edges_assignment: Some(ti_vec![Boundary, Valley, Boundary, Boundary, Boundary, Boundary, Boundary]),
            edges_fold_angle_f64: Some(ti_vec![0.0, -180.0, 0.0, 0.0, 0.0, 0.0, 0.0]),
            edges_fold_angle_exact: Some(ti_vec![
                Angle::ZERO, -Angle::DEG_180, Angle::ZERO, Angle::ZERO, Angle::ZERO, Angle::ZERO, Angle::ZERO,
            ]),
            edges_length_f64: Some(ti_vec![1.0, 1.0, 1.0, 2.0, 2.0, 2.0, 2.0]),
            edges_length2_exact: Some(ti_vec![
                based_expr!(1), based_expr!(1), based_expr!(1), based_expr!(2), based_expr!(2), based_expr!(2), based_expr!(2)
            ]),
            edge_orders: Some(vec![
                (E(0), E(1), EdgeOrder::Right),
                (E(1), E(0), EdgeOrder::Left),
                (E(1), E(2), EdgeOrder::Right),
                (E(2), E(1), EdgeOrder::Left),
            ]),
            edges_custom: indexmap! {
                "test:test".to_owned() => ti_vec![json!(1), json!(2), json!(3), json!(4), json!(5), json!(6), json!(7)]
            },
            ..Default::default()
        }.with_topology_vh_fh(Some(ti_vec![
            vec![H(0), H(6)],
            vec![H(7), H(2), H(8)],
            vec![H(9), H(4)],
            vec![H(10), H(1)],
            vec![H(12), H(3), H(11)],
            vec![H(5), H(13)],
        ]), Some(ti_vec![
            vec![H(0), H(10), H(3), H(7)],
            vec![H(2), H(12), H(5), H(9)],
        ]));

        let mut frame = wrap_fold.clone();
        let (edge_data, face_datas) = frame.swap_remove_edge(E(0));
        assert_eq!(edge_data, EdgeData {
            vertices: [V(0), V(3)],
            face_corners: Some([vec![], vec![]]),
            assignment: Some(EdgeAssignment::Boundary),
            fold_angle_f64: Some(0.0),
            fold_angle_exact: Some(Angle::ZERO),
            length_f64: Some(1.0),
            length2_exact: Some(based_expr!(1)),
            custom: indexmap! { "test:test".to_owned() => json!(1) }
        });
        assert_eq!(face_datas, indexmap! {
            F(0) => FaceData {
                half_edges: vec![H(0), H(10), H(3), H(7)],
                custom: indexmap! {}
            }
        });
        assert_eq!(frame.vertices_half_edges, Some(ti_vec![
            vec![H(6)],
            vec![H(7), H(2), H(8)],
            vec![H(9), H(4)],
            vec![H(10)],
            vec![H(0), H(3), H(11)],
            vec![H(5), H(1)],
        ]));
        assert_eq!(frame.faces_half_edges, Some(ti_vec![
            vec![H(2), H(0), H(5), H(9)],
        ]));
        assert_eq!(frame.edges_vertices, Some(ti_vec![
            [V(4), V(5)],
            [V(1), V(4)],
            [V(2), V(5)],
            [V(0), V(1)],
            [V(1), V(2)],
            [V(3), V(4)],
        ]));
        assert_eq!(frame.edges_face_corners, Some(ti_vec![
            [vec![C(F(0), 1)], vec![]],
            [vec![C(F(0), 0)], vec![]],
            [vec![], vec![C(F(0), 2)]],
            [vec![], vec![]],
            [vec![], vec![C(F(0), 3)]],
            [vec![], vec![]],
        ]));
        assert_eq!(frame.edges_assignment, Some(ti_vec![Boundary, Valley, Boundary, Boundary, Boundary, Boundary]));
        assert_eq!(frame.edges_fold_angle_f64, Some(ti_vec![0.0, -180.0, 0.0, 0.0, 0.0, 0.0]));
        assert_eq!(frame.edges_fold_angle_exact, Some(ti_vec![
            Angle::ZERO, -Angle::DEG_180, Angle::ZERO, Angle::ZERO, Angle::ZERO, Angle::ZERO
        ]));
        assert_eq!(frame.edges_length_f64, Some(ti_vec![2.0, 1.0, 1.0, 2.0, 2.0, 2.0]));
        assert_eq!(frame.edges_length2_exact, Some(ti_vec![
            based_expr!(2), based_expr!(1), based_expr!(1), based_expr!(2), based_expr!(2), based_expr!(2)
        ]));
        assert_eq!(frame.edge_orders, Some(vec![
            (E(1), E(2), EdgeOrder::Right),
            (E(2), E(1), EdgeOrder::Left),
        ]));
        assert_eq!(frame.edges_custom, indexmap! {
            "test:test".to_owned() => ti_vec![json!(7), json!(2), json!(3), json!(4), json!(5), json!(6)]
        });

        // Remove the edge adjacent to both faces
        let mut frame = wrap_fold.clone();
        let (edge_data, face_datas) = frame.swap_remove_edge(E(1));
        assert_eq!(edge_data, EdgeData {
            vertices: [V(1), V(4)],
            face_corners: Some([vec![], vec![]]),
            assignment: Some(EdgeAssignment::Valley),
            fold_angle_f64: Some(-180.0),
            fold_angle_exact: Some(-Angle::DEG_180),
            length_f64: Some(1.0),
            length2_exact: Some(based_expr!(1)),
            custom: indexmap! { "test:test".to_owned() => json!(2) }
        });
        assert_eq!(face_datas, indexmap! {
            F(0) => FaceData {
                half_edges: vec![H(0), H(10), H(3), H(7)],
                custom: indexmap! {}
            },
            F(1) => FaceData {
                half_edges: vec![H(2), H(12), H(5), H(9)],
                custom: indexmap! {}
            }
        });
        assert_eq!(frame.vertices_half_edges, Some(ti_vec![
            vec![H(0), H(6)],
            vec![H(7), H(8)],
            vec![H(9), H(4)],
            vec![H(10), H(1)],
            vec![H(2),  H(11)],
            vec![H(5), H(3)],
        ]));
        assert_eq!(frame.faces_half_edges, Some(ti_vec![]));
        assert_eq!(frame.edges_vertices, Some(ti_vec![
            [V(0), V(3)],
            [V(4), V(5)],
            [V(2), V(5)],
            [V(0), V(1)],
            [V(1), V(2)],
            [V(3), V(4)],
        ]));
        assert_eq!(frame.edges_face_corners, Some(ti_vec![
            [vec![], vec![]],
            [vec![], vec![]],
            [vec![], vec![]],
            [vec![], vec![]],
            [vec![], vec![]],
            [vec![], vec![]],
        ]));
        assert_eq!(frame.edges_assignment, Some(ti_vec![Boundary, Boundary, Boundary, Boundary, Boundary, Boundary]));
        assert_eq!(frame.edges_fold_angle_f64, Some(ti_vec![0.0, 0.0, 0.0, 0.0, 0.0, 0.0]));
        assert_eq!(frame.edges_fold_angle_exact, Some(ti_vec![
            Angle::ZERO, Angle::ZERO, Angle::ZERO, Angle::ZERO, Angle::ZERO, Angle::ZERO
        ]));
        assert_eq!(frame.edges_length_f64, Some(ti_vec![1.0, 2.0, 1.0, 2.0, 2.0, 2.0]));
        assert_eq!(frame.edges_length2_exact, Some(ti_vec![
            based_expr!(1), based_expr!(2), based_expr!(1), based_expr!(2), based_expr!(2), based_expr!(2)
        ]));
        assert_eq!(frame.edge_orders, Some(vec![]));
        assert_eq!(frame.edges_custom, indexmap! {
            "test:test".to_owned() => ti_vec![json!(1), json!(7), json!(3), json!(4), json!(5), json!(6)]
        });

        // Slitted face, remove edge adjacent to the same face multiple times
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
        let (edge_data, face_datas) = frame.swap_remove_edge(E(8));
        assert_eq!(edge_data, EdgeData {
            vertices: [V(0), V(4)],
            face_corners: Some([vec![], vec![]]),
            assignment: None,
            fold_angle_f64: None,
            fold_angle_exact: None,
            length_f64: None,
            length2_exact: None,
            custom: indexmap! {}
        });
        assert_eq!(face_datas, indexmap! {
            F(0) => FaceData {
                half_edges: vec![H(0), H(2), H(4), H(6), H(16), H(15), H(13), H(11), H(9), H(17)],
                custom: indexmap! {}
            },
        });
        assert_eq!(frame.vertices_half_edges, Some(ti_vec![
            vec![H(0), H(7)],
            vec![H(1), H(2)],
            vec![H(3), H(4)],
            vec![H(5), H(6)],
            vec![H(8), H(15)],
            vec![H(9), H(10)],
            vec![H(11), H(12)],
            vec![H(13), H(14)],
        ]));
        assert_eq!(frame.faces_half_edges, Some(ti_vec![
            vec![H(8), H(10), H(12), H(14)],
        ]));
        assert_eq!(frame.edges_vertices, Some(ti_vec![
            [V(0), V(1)],
            [V(1), V(2)],
            [V(2), V(3)],
            [V(3), V(0)],
            [V(4), V(5)],
            [V(5), V(6)],
            [V(6), V(7)],
            [V(7), V(4)],
        ]));
        assert_eq!(frame.edges_face_corners, Some(ti_vec![
            [vec![], vec![]],
            [vec![], vec![]],
            [vec![], vec![]],
            [vec![], vec![]],
            [vec![C(F(0), 0)], vec![]],
            [vec![C(F(0), 1)], vec![]],
            [vec![C(F(0), 2)], vec![]],
            [vec![C(F(0), 3)], vec![]],
        ]));
    }

    #[test]
    fn test_shift_remove_face() {
        // 0---4---1---5---2---6---3
        // |       |       |       |
        // 0   0   1   1   2   2   3
        // |       |       |       |
        // 4---7---5---8---6---9---7
        // all edges point down and right
        // all faces start at the top-left corner and go ccw
        let wrap_fold = Frame {
            face_orders: Some(vec![
                (F(0), F(1), FaceOrder::Above),
                (F(1), F(0), FaceOrder::Above),
                (F(0), F(2), FaceOrder::Above),
                (F(2), F(0), FaceOrder::Below),
                (F(2), F(1), FaceOrder::Above),
                (F(1), F(2), FaceOrder::Above),
            ]),
            faces_custom: indexmap! {
                "test:test".to_owned() => ti_vec![json!(1), json!(2), json!(4)]
            },
            ..Default::default()
        }.with_topology_vh_fh(Some(ti_vec![
            vec![H(0), H(8)],
            vec![H(9), H(2), H(10)],
            vec![H(11), H(4), H(12)],
            vec![H(13), H(6)],
            vec![H(14), H(1)],
            vec![H(16), H(3), H(15)],
            vec![H(18), H(5), H(17)],
            vec![H(7), H(19)],
        ]), Some(ti_vec![
            vec![H(0), H(14), H(3), H(9)],
            vec![H(2), H(16), H(5), H(11)],
            vec![H(4), H(18), H(7), H(13)],
        ]));

        let mut frame = wrap_fold.clone();
        let face_data = frame.shift_remove_face(F(0));
        assert_eq!(face_data, FaceData {
            half_edges: vec![H(0), H(14), H(3), H(9)],
            custom: indexmap! { "test:test".to_owned() => json!(1) }
        });
        assert_eq!(frame.faces_half_edges, Some(ti_vec![
            vec![H(2), H(16), H(5), H(11)],
            vec![H(4), H(18), H(7), H(13)],
        ]));
        assert_eq!(frame.edges_face_corners, Some(ti_vec![
            [vec![], vec![]],
            [vec![C(F(0), 0)], vec![]],
            [vec![C(F(1), 0)], vec![C(F(0), 2)]],
            [vec![], vec![C(F(1), 2)]],
            [vec![], vec![]],
            [vec![], vec![C(F(0), 3)]],
            [vec![], vec![C(F(1), 3)]],
            [vec![], vec![]],
            [vec![C(F(0), 1)], vec![]],
            [vec![C(F(1), 1)], vec![]],
        ]));
        assert_eq!(frame.face_orders, Some(vec![
            (F(1), F(0), FaceOrder::Above),
            (F(0), F(1), FaceOrder::Above),
        ]));
        assert_eq!(frame.faces_custom, indexmap! {
            "test:test".to_owned() => ti_vec![json!(2), json!(4)]
        });

        // Removing the last entry
        let mut frame = wrap_fold.clone();
        let face_data = frame.shift_remove_face(F(2));
        assert_eq!(face_data, FaceData {
            half_edges: vec![H(4), H(18), H(7), H(13)],
            custom: indexmap! { "test:test".to_owned() => json!(4) }
        });
        assert_eq!(frame.faces_half_edges, Some(ti_vec![
            vec![H(0), H(14), H(3), H(9)],
            vec![H(2), H(16), H(5), H(11)],
        ]));
        assert_eq!(frame.edges_face_corners, Some(ti_vec![
            [vec![C(F(0), 0)], vec![]],
            [vec![C(F(1), 0)], vec![C(F(0), 2)]],
            [vec![], vec![C(F(1), 2)]],
            [vec![], vec![]],
            [vec![], vec![C(F(0), 3)]],
            [vec![], vec![C(F(1), 3)]],
            [vec![], vec![]],
            [vec![C(F(0), 1)], vec![]],
            [vec![C(F(1), 1)], vec![]],
            [vec![], vec![]],
        ]));
        assert_eq!(frame.face_orders, Some(vec![
            (F(0), F(1), FaceOrder::Above),
            (F(1), F(0), FaceOrder::Above),
        ]));
        assert_eq!(frame.faces_custom, indexmap! {
            "test:test".to_owned() => ti_vec![json!(1), json!(2)]
        });
    }

    #[test]
    fn test_swap_remove_face() {
        // 0---4---1---5---2---6---3
        // |       |       |       |
        // 0   0   1   1   2   2   3
        // |       |       |       |
        // 4---7---5---8---6---9---7
        // all edges point down and right
        // all faces start at the top-left corner and go ccw
        let wrap_fold = Frame {
            face_orders: Some(vec![
                (F(0), F(1), FaceOrder::Above),
                (F(1), F(0), FaceOrder::Above),
                (F(0), F(2), FaceOrder::Above),
                (F(2), F(0), FaceOrder::Below),
                (F(2), F(1), FaceOrder::Above),
                (F(1), F(2), FaceOrder::Above),
            ]),
            faces_custom: indexmap! {
                "test:test".to_owned() => ti_vec![json!(1), json!(2), json!(4)]
            },
            ..Default::default()
        }.with_topology_vh_fh(Some(ti_vec![
            vec![H(0), H(8)],
            vec![H(9), H(2), H(10)],
            vec![H(11), H(4), H(12)],
            vec![H(13), H(6)],
            vec![H(14), H(1)],
            vec![H(16), H(3), H(15)],
            vec![H(18), H(5), H(17)],
            vec![H(7), H(19)],
        ]), Some(ti_vec![
            vec![H(0), H(14), H(3), H(9)],
            vec![H(2), H(16), H(5), H(11)],
            vec![H(4), H(18), H(7), H(13)],
        ]));

        let mut frame = wrap_fold.clone();
        let face_data = frame.swap_remove_face(F(0));
        assert_eq!(face_data, FaceData {
            half_edges: vec![H(0), H(14), H(3), H(9)],
            custom: indexmap! { "test:test".to_owned() => json!(1) }
        });
        assert_eq!(frame.faces_half_edges, Some(ti_vec![
            vec![H(4), H(18), H(7), H(13)],
            vec![H(2), H(16), H(5), H(11)],
        ]));
        assert_eq!(frame.edges_face_corners, Some(ti_vec![
            [vec![], vec![]],
            [vec![C(F(1), 0)], vec![]],
            [vec![C(F(0), 0)], vec![C(F(1), 2)]],
            [vec![], vec![C(F(0), 2)]],
            [vec![], vec![]],
            [vec![], vec![C(F(1), 3)]],
            [vec![], vec![C(F(0), 3)]],
            [vec![], vec![]],
            [vec![C(F(1), 1)], vec![]],
            [vec![C(F(0), 1)], vec![]],
        ]));
        assert_eq!(frame.face_orders, Some(vec![
            (F(0), F(1), FaceOrder::Above),
            (F(1), F(0), FaceOrder::Above),
        ]));
        assert_eq!(frame.faces_custom, indexmap! {
            "test:test".to_owned() => ti_vec![json!(4), json!(2)]
        });

        // Removing the last entry
        let mut frame = wrap_fold.clone();
        let face_data = frame.swap_remove_face(F(2));
        assert_eq!(face_data, FaceData {
            half_edges: vec![H(4), H(18), H(7), H(13)],
            custom: indexmap! { "test:test".to_owned() => json!(4) }
        });
        assert_eq!(frame.faces_half_edges, Some(ti_vec![
            vec![H(0), H(14), H(3), H(9)],
            vec![H(2), H(16), H(5), H(11)],
        ]));
        assert_eq!(frame.edges_face_corners, Some(ti_vec![
            [vec![C(F(0), 0)], vec![]],
            [vec![C(F(1), 0)], vec![C(F(0), 2)]],
            [vec![], vec![C(F(1), 2)]],
            [vec![], vec![]],
            [vec![], vec![C(F(0), 3)]],
            [vec![], vec![C(F(1), 3)]],
            [vec![], vec![]],
            [vec![C(F(0), 1)], vec![]],
            [vec![C(F(1), 1)], vec![]],
            [vec![], vec![]],
        ]));
        assert_eq!(frame.face_orders, Some(vec![
            (F(0), F(1), FaceOrder::Above),
            (F(1), F(0), FaceOrder::Above),
        ]));
        assert_eq!(frame.faces_custom, indexmap! {
            "test:test".to_owned() => ti_vec![json!(1), json!(2)]
        });
    }

    #[test]
    fn test_flip_face() {
        // 0---3---1---4---2
        // |       |       |
        // 0   0   1   1   2
        // |       |       |
        // 3---5---4---6---5
        // all edges point down and right
        // all faces start at the top-left corner and go ccw
        let mut frame = Frame {
            ..Default::default()
        }.with_topology_vh_fh(Some(ti_vec![
            vec![H(0), H(6)],
            vec![H(7), H(2), H(8)],
            vec![H(9), H(4)],
            vec![H(10), H(1)],
            vec![H(12), H(3), H(11)],
            vec![H(5), H(13)],
        ]), Some(ti_vec![
            vec![H(0), H(10), H(3), H(7)],
            vec![H(2), H(12), H(5), H(9)],
        ]));

        frame.flip_face(F(0));
        frame.assert_topology_consistent();
        assert_eq!(frame.faces_half_edges, Some(ti_vec![
            vec![H(6), H(2), H(11), H(1)],
            vec![H(2), H(12), H(5), H(9)],
        ]));

        frame.flip_face(F(1));
        frame.assert_topology_consistent();
        assert_eq!(frame.faces_half_edges, Some(ti_vec![
            vec![H(6), H(2), H(11), H(1)],
            vec![H(8), H(4), H(13), H(3)],
        ]));

        // Slitted face
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

        frame.flip_face(F(0));
        frame.assert_topology_consistent();
        assert_eq!(frame.faces_half_edges, Some(ti_vec![
            vec![H(16), H(8), H(10), H(12), H(14), H(17), H(7), H(5), H(3), H(1)],
            vec![H(8), H(10), H(12), H(14)],
        ]));

        // A möbius face. Noticably, a face duplicates a half-edge.
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

        frame.flip_face(F(0));
        frame.assert_topology_consistent();
        assert_eq!(frame.faces_half_edges, Some(ti_vec![
            vec![H(8), H(4), H(6), H(8), H(3), H(1)],
            vec![H(0), H(2), H(4), H(6)],
        ]));
    }
}