use std::{fmt::{Debug, Display}, fs::File, io::BufReader, ops::Index, path::Path};
use derive_more::{Add, AddAssign, Div, DivAssign, From, Into, Mul, MulAssign, Rem, RemAssign, Sub, SubAssign};
use exact_number::{basis::ArcBasis, BasedExpr};
use indexmap::IndexMap;
use nalgebra::DMatrix;
use serde::{Deserialize, Serialize};
use serde_json::{Map, Value};
use serde_repr::{Deserialize_repr, Serialize_repr};
use typed_index_collections::{TiSlice, TiVec};

/// A subjective interpretation about what the entire file represents.
#[derive(Clone, Debug, PartialEq, Eq, PartialOrd, Ord)]
#[derive(Serialize, Deserialize)]
#[serde(rename_all = "camelCase")]
pub enum FileClass {
    /// A single origami model, possibly still in multiple frames to represent crease pattern, folded form, etc.
    SingleModel,
    /// Multiple origami models collected together into one file
    MultiModel,
    /// Animation of sequence of frames, e.g., illustrating a continuous folding motion
    Animation,
    /// A sequence of frames representing folding steps, as in origami diagrams
    Diagrams,
    #[serde(untagged)]
    Custom(String),
}

pub(crate) fn version() -> String { "1.2".to_owned() }

#[derive(Clone, Debug)]
#[cfg_attr(test, derive(PartialEq))] // Really no point in this derivation outside of tests
pub struct Fold {
    /// The version of the FOLD spec that the file assumes.
    /// **Strongly recommended**, in case we ever have to make
    /// backward-incompatible changes.
    pub file_spec: String,
    /// The software that created the file.
    /// **Recommended** for files output by computer software;
    /// less important for files made by hand.
    pub file_creator: Option<String>,
    /// The human author
    pub file_author: Option<String>,
    /// A title for the entire file.
    pub file_title: Option<String>,
    /// A description of the entire file.
    pub file_description: Option<String>,
    /// A subjective interpretation about what the entire file represents.
    pub file_classes: Vec<FileClass>,
    /// The frames in the file. The key frame is frame 0.
    pub file_frames: TiVec<FrameIndex, Frame>,
    /// Custom data
    pub file_custom: IndexMap<String, Value>,
}

impl Fold {
    /// Gets the key frame, a.k.a. frame 0.
    pub fn key_frame(&self) -> &Frame {
        self.file_frames.first().unwrap()
    }

    /// Gets the key frame, a.k.a. frame 0, mutably
    pub fn key_frame_mut(&mut self) -> &mut Frame {
        self.file_frames.first_mut().unwrap()
    }
}

impl Default for Fold {
    fn default() -> Self {
        Self {
            file_spec: version(),
            file_creator: Default::default(),
            file_author: Default::default(),
            file_title: Default::default(),
            file_description: Default::default(),
            file_classes: Default::default(),
            file_frames: Default::default(),
            file_custom: Default::default(),
        }
    }
}

/// A subjective interpretation about what the frame represents.
#[derive(Clone, Debug, PartialEq, Eq, PartialOrd, Ord)]
#[derive(Serialize, Deserialize)]
#[serde(rename_all = "camelCase")]
pub enum FrameClass {
    /// a crease pattern (unfolded)
    CreasePattern,
    /// a folded form/state, e.g. flat folding or 3D folding
    FoldedForm,
    /// vertices and edges, but no lengths or faces
    Graph,
    /// vertices and edges and edge lengths, but no faces
    Linkage,
    #[serde(untagged)]
    Custom(String),
}

/// Attribute that objectively describes a property of the
/// folded structure being represented.
#[derive(Clone, Debug, PartialEq, Eq, PartialOrd, Ord)]
#[derive(Serialize, Deserialize)]
#[serde(rename_all = "camelCase")]
pub enum FrameAttribute {
    /// the coordinates lie in 2D (xy); z coordinates are all implicitly or explicitly 0
    #[serde(rename = "2D")]
    _2D,
    /// the coordinates lie in 3D (xyz) and not 2D (xy)
    #[serde(rename = "3D")]
    _3D,
    /// the polyhedral complex is not embedded in Euclidean space,
    /// so there are no vertex coordinates (but there might be edge lengths defining intrinsic geometry)
    Abstract,
    /// the polyhedral complex is a manifold (has at most two faces incident to each edge)
    Manifold,
    /// the polyhedral complex is not a manifold
    NonManifold,
    /// the polyhedral complex is orientable, meaning it can be assigned a consistent normal direction (and hence it is also manifold)
    Orientable,
    /// the polyhedral complex is not orientable, meaning it cannot be assigned a consistent normal direction
    NonOrientable,
    /// the polyhedral complex has faces that touch in their relative interiors, in which case you probably want a face ordering
    SelfTouching,
    /// the polyhedral complex does not have faces that touch in their relative interiors
    NonSelfTouching,
    /// the polyhedral complex has properly intersecting faces
    SelfIntersecting,
    /// the polyhedral complex does not have properly intersecting faces
    NonSelfIntersecting,
    /// an edge has an assignment of `Cut` (cut/slit representing multiple `Boundary` edges)
    Cuts,
    /// no edges have an assignment of `Cut` (cut/slit representing multiple `Boundary` edges)
    NoCuts,
    /// an edge has an assignment of `Join` (join)
    Joins,
    /// no edges have an assignment of `Join` (join)
    NoJoins,
    /// all faces are convex polygons
    ConvexFaces,
    /// some faces are nonconvex
    NonConvexFaces,
    #[serde(untagged)]
    Custom(String),
}

/// Physical or logical unit that all coordinates are relative to.
#[derive(Clone, Debug, PartialEq, Eq, PartialOrd, Ord)]
#[derive(Serialize, Deserialize)]
pub enum FrameUnit {
    /// (equivalent to not specifying a unit): no physical meaning
    #[serde(rename = "unit")]
    Unit,
    /// inches (25.4 mm)
    #[serde(rename = "in")]
    Inch,
    /// desktop publishing/PostScript [points](https://en.wikipedia.org/wiki/Point_(typography)) (1/72 in)
    #[serde(rename = "pt")]
    Point,
    /// meters (1/299,792,458 light seconds)
    #[serde(rename = "m")]
    Meter,
    /// centimeters (1/100 meters)
    #[serde(rename = "cm")]
    Centimeter,
    /// millimeters (1/1000 meters)
    #[serde(rename = "mm")]
    Millimeter,
    /// microns (1/1,000,000 meters)
    #[serde(rename = "um")]
    Micrometer,
    /// nanometers (1/1,000,000,000 meters)
    #[serde(rename = "nm")]
    Nanometer,
    /// Beware, this will probably not be understood by software
    #[serde(untagged)]
    Custom(String)
}

impl Default for FrameUnit {
    fn default() -> Self {
        Self::Unit
    }
}

impl FrameUnit {
    pub fn is_default(&self) -> bool { self == &Self::Unit }
}

/// For each edge, a string representing its fold direction assignment
#[derive(Clone, Copy, Debug, PartialEq, Eq, PartialOrd, Ord)]
#[derive(Serialize, Deserialize)]
pub enum EdgeAssignment {
    /// border/boundary edge (only one incident face)
    #[serde(rename = "B")]
    Boundary,
    /// mountain crease
    #[serde(rename = "M")]
    Mountain,
    /// valley crease
    #[serde(rename = "V")]
    Valley,
    /// flat (unfolded) crease
    #[serde(rename = "F")]
    Flat,
    /// unassigned/unknown crease
    #[serde(rename = "U")]
    Unassigned,
    /// cut/slit edge (should be treated as multiple `Boundary` edges)
    #[serde(rename = "C")]
    Cut,
    /// join edge (incident faces should be treated as a single face)
    #[serde(rename = "J")]
    Join,
}

#[derive(Clone, Copy, Debug, PartialEq, Eq, PartialOrd, Ord)]
#[derive(Serialize_repr, Deserialize_repr)]
#[repr(i8)]
pub enum FaceOrder {
    /// Given `(f, g, Below)`, face `f` lies *below face `g`, i.e., on the side opposite `g`'s normal vector in the folded state.
    Below = -1,
    /// Given `(f, g, Unknown)`, `f` and `g` have unknown stacking order (e.g., they do not overlap in their interiors).
    Unknown = 0,
    /// Given `(f, g, Above)`, face `f` lies *above* face `g`, i.e., on the side pointed to by `g`'s normal vector in the folded state.
    Above = 1,
}

#[derive(Clone, Copy, Debug, PartialEq, Eq, PartialOrd, Ord)]
#[derive(Serialize_repr, Deserialize_repr)]
#[repr(i8)]
pub enum EdgeOrder {
    /// Given `(e, f, Right)`, edge `e` lies locally on the *right* side of edge `f` (relative to edge `f`'s orientation given by `edges_vertices`)
    Right = -1,
    /// Given `(e, f, Unknown)`, `e` and `f` have unknown stacking order (e.g., they do not overlap in their interiors).
    Unknown = 0,
    /// Given `(e, f, Left)`, edge `e` lies locally on the *left* side of edge `f` (relative to edge `f`'s orientation given by `edges_vertices`)
    Left = 1,
}

/// A matrix of coordinates; either exact or approximate.
#[derive(Clone, Copy, Debug, PartialEq)]
pub enum CoordsRef<'a> {
    Exact(&'a DMatrix<BasedExpr>),
    Approx(&'a DMatrix<f64>)
}

impl<'a> CoordsRef<'a> {
    /// Gets the number of vertices.
    pub fn num_vertices(&self) -> usize {
        match self {
            CoordsRef::Exact(c) => c.ncols(),
            CoordsRef::Approx(c) => c.ncols(),
        }
    }

    /// Gets the number of dimensions.
    pub fn num_dimensions(&self) -> usize {
        match self {
            CoordsRef::Exact(c) => c.nrows(),
            CoordsRef::Approx(c) => c.nrows(),
        }
    }
}

macro_rules! impl_display_index {
    (impl Debug, Display for $ty:ty { $fmt:literal }) => {
        impl std::fmt::Debug for $ty {
            fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
                write!(f, $fmt, self.0)
            }
        }

        impl std::fmt::Display for $ty {
            fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
                <$ty as Debug>::fmt(&self, f)
            }
        }
    };
}

/// Type-safe vertex index.
/// However, for coordinates, you need to get the `usize`.
#[derive(Clone, Copy, PartialEq, Eq, Hash, PartialOrd, Ord, From, Into,
    Add, AddAssign, Sub, SubAssign, Mul, MulAssign, Div, DivAssign, Rem, RemAssign)]
#[derive(Serialize, Deserialize)]
#[repr(transparent)]
#[serde(transparent)]
pub struct Vertex(pub usize);
impl_display_index! { impl Debug, Display for Vertex { "v_{}" } }

/// Type-safe edge index.
#[derive(Clone, Copy, PartialEq, Eq, Hash, PartialOrd, Ord, From, Into,
    Add, AddAssign, Sub, SubAssign, Mul, MulAssign, Div, DivAssign, Rem, RemAssign)]
#[derive(Serialize, Deserialize)]
#[repr(transparent)]
#[serde(transparent)]
pub struct Edge(pub usize);
impl_display_index! { impl Debug, Display for Edge { "e_{}" } }

/// Type-safe half-edge index.
/// 
/// A *half-edge* is an edge with a direction.
/// The edge index is the half-edge index divided by 2,
/// and the direction is the direction specified by `edges_vertices` if even
/// and the opposite direction if `odd`.
#[derive(Clone, Copy, PartialEq, Eq, Hash, PartialOrd, Ord, From, Into,
    Add, AddAssign, Sub, SubAssign, Mul, MulAssign, Div, DivAssign, Rem, RemAssign)]
#[derive(Serialize, Deserialize)]
#[repr(transparent)]
#[serde(transparent)]
pub struct HalfEdge(pub usize);

impl HalfEdge {
    /// Constructs a half-edge index with an edge index
    pub fn new(edge: Edge, flipped: bool) -> Self {
        Self(edge.0 * 2 + (flipped as usize))
    }

    /// Get the value of the flip bit.
    pub fn flip_bit(self) -> bool {
        (self.0 & 1) != 0
    }

    /// Get the edge.
    pub fn edge(self) -> Edge {
        Edge(self.0 >> 1)
    }

    /// Flip the half-edge.
    pub fn flip(&mut self) {
        self.0 ^= 1;
    }

    /// Get the half-edge flipped.
    pub fn flipped(mut self) -> Self {
        self.flip();
        self
    }
}

impl Debug for HalfEdge {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        write!(f, "h_{}{}", if self.0 % 2 == 0 {">"} else {"<"}, self.0 / 2)
    }
}

impl Display for HalfEdge {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        <HalfEdge as Debug>::fmt(&self, f)
    }
}

/// Type-safe face index.
#[derive(Clone, Copy, PartialEq, Eq, Hash, PartialOrd, Ord, From, Into,
    Add, AddAssign, Sub, SubAssign, Mul, MulAssign, Div, DivAssign, Rem, RemAssign)]
#[derive(Serialize, Deserialize)]
#[repr(transparent)]
#[serde(transparent)]
pub struct Face(pub usize);
impl_display_index! { impl Debug, Display for Face { "f_{}" } }

/// Type-safe frame index.
#[derive(Clone, Copy, Debug, PartialEq, Eq, Hash, PartialOrd, Ord, From, Into,
    Add, AddAssign, Sub, SubAssign, Mul, MulAssign, Div, DivAssign, Rem, RemAssign)]
#[derive(Serialize, Deserialize)]
#[repr(transparent)]
#[serde(transparent)]
pub struct FrameIndex(pub usize);

trait AtHalfEdge {
    type Output;
    // because `std::ops::Index` forces a reference to be returned
    fn at(&self, index: HalfEdge) -> Self::Output;
}

pub type EdgeVerticesSlice = TiSlice<Edge, [Vertex; 2]>;

impl AtHalfEdge for EdgeVerticesSlice {
    type Output = [Vertex; 2];

    fn at(&self, index: HalfEdge) -> Self::Output {
        let mut result = self[index.edge()];
        if index.flip_bit() { result.reverse() }
        result
    }
}

/// A FOLD frame, containing geometry information.
/// 
/// All operations prefer exact coordinates (`vertices_coords_exact`, `edges_fold_angle_cs`, `edges_length2_exact`),
/// not messing with approximate coordinates (`vertices_coords_f64`, `edges_fold_angle_f64`, `edges_length_f64`) when both exist.
#[derive(Clone, Debug)]
#[derive(Serialize, Deserialize)]
#[cfg_attr(test, derive(PartialEq))] // Really no point in this derivation outside of tests
#[serde(try_from = "crate::ser_de::SerDeFrame", into = "crate::ser_de::SerDeFrame")]
pub struct Frame {
    /// The human author
    pub frame_author: Option<String>,
    /// A title for the frame.
    pub frame_title: Option<String>,
    /// A description of the frame.
    pub frame_description: Option<String>,
    /// A subjective interpretation about what the frame represents.
    pub frame_classes: Vec<FrameClass>,
    /// Attributes that objectively describe properties of the
    /// folded structure being represented.
    /// 
    /// # Requirements
    /// * At most 1 of `_2D`, `_3D`, and `_Abstract` are set.
    /// * `Manifold` and `NotManifold` are not both set.
    /// * `Orientable` and `NotManifold` are not both set.
    /// * `Orientable` and `NotOrientable` are not both set.
    /// * `SelfTouching` and `NotSelfTouching` are not both set.
    /// * `SelfIntersecting` and `NotSelfIntersecting` are not both set.
    /// * `Cuts` and `NoCuts` are not both set.
    /// * `Joins` and `NoJoints` are not both set.
    /// * `ConvexFaces` and `NoConvexFaces` are not both set.
    pub frame_attributes: Vec<FrameAttribute>,
    /// Physical or logical unit that all coordinates are relative to
    pub frame_unit: FrameUnit,

    /// The basis that exact coordinates uses.
    pub basis: Option<ArcBasis>,

    /// For each vertex, an array of approximate coordinates,
    /// such as `[x, y, z]` or `[x, y]` (where `z` is implicitly zero).
    /// In higher dimensions, all trailing unspecified coordinates are implicitly
    /// zero.
    /// Vertex coordinates are columns.
    /// 
    /// **Recommended** except for frames with attribute `Abstract`.
    /// 
    /// # Requirements
    /// * **Exists** if `FrameAttribute::_2D` or `FrameAttribute::_3D` is set
    ///   and `vertices_coords_exact == None`.
    /// * If `FrameAttribute::_2D` or `FrameAttribute::_3D` is set,
    ///   the coordinate dimensions match the attribute.
    pub vertices_coords_f64: Option<DMatrix<f64>>,
    /// For each vertex, an array of exact coordinates,
    /// such as `[x, y, z]` or `[x, y]` (where `z` is implicitly zero).
    /// In higher dimensions, all trailing unspecified coordinates are implicitly
    /// zero.
    /// Vertex coordinates are columns.
    /// 
    /// **Recommended** except for frames with attribute `Abstract`.
    /// 
    /// # Requirements
    /// * If `FrameAttribute::_2D` or `FrameAttribute::_3D` is set,
    ///   the coordinate dimensions match the attribute.
    pub vertices_coords_exact: Option<DMatrix<BasedExpr>>,
    /// The number of vertices. Necessary to handle isolated vertices.
    pub num_vertices: usize,
    /// For each vertex, an array of edge IDs for the *half-edge*s
    /// incident to the vertex.  If the frame represents an orientable manifold,
    /// this list should be ordered counterclockwise around the vertex.
    /// If the frame is a nonorientable manifold, this list should be cyclically
    /// ordered around the vertex.
    /// Note that `vertices_half_edges[v][i]` should be an edge *from* vertex
    /// `v` *to* the edge's other vertex.
    /// 
    /// This should be calculated from `edges_vertices`.
    pub vertices_half_edges: Option<TiVec<Vertex, Vec<HalfEdge>>>,
    /// `edges_vertices`: For each edge, an array `[u, v]` of two vertex IDs for
    /// the two endpoints of the edge.  This effectively defines the *orientation*
    /// of the edge, from `u` to `v`.  (This orientation choice is arbitrary,
    /// but is used to define the ordering of `edges_faces`.)
    /// **Recommended** in frames having any `edges_...` property
    /// (e.g., to represent mountain-valley assignment).
    /// 
    /// # Requirements
    /// * No edge can repeat a vertex. (no self-loops)
    pub edges_vertices: Option<TiVec<Edge, [Vertex; 2]>>,
    /// For each edge, an array of face IDs for the faces incident
    /// to the edge, possibly including `None`s.
    /// For nonmanifolds in particular, the `Some()` faces should be listed in
    /// counterclockwise order around the edge,
    /// relative to the orientation of the edge.
    /// For manifolds, the array for each edge should be an array of length 2,
    /// where the first entry is the face locally to the "left" of the edge
    /// (or `None` if there is no such face) and the second entry is the face
    /// locally to the "right" of the edge (or `None` if there is no such face);
    /// for orientable manifolds, "left" and "right" must be consistent with the
    /// manifold orientation given by the counterclockwise orientation of faces.
    /// However, a boundary edge may also be represented by a length-1 array, with
    /// the `None` omitted, to be consistent with the nonmanifold representation.
    /// 
    /// This should be calculated from `faces_edges`.
    pub edges_faces: Option<TiVec<Edge, Vec<Option<Face>>>>,
    /// For each edge, a string representing its fold direction assignment
    pub edges_assignment: Option<TiVec<Edge, EdgeAssignment>>,
    /// For each edge, the fold angle (deviation from flatness)
    /// along each edge of the pattern.  The fold angle is a number in degrees
    /// lying in the range [-180, 180].  The fold angle is positive for
    /// valley folds, negative for mountain folds, and zero for flat, unassigned,
    /// and border folds.  Accordingly, the sign of `edge_foldAngle` should match
    /// `edges_assignment` if both are specified.
    /// 
    /// For now, this is an `f64` regardless of `N`, because finding a good representation
    /// with `N` is tricky. `[cos(a), sin(a)]` doesn't work because -180° and 180° are both in range.
    pub edges_fold_angle_f64: Option<TiVec<Edge, f64>>,
    /// For each edge, the *exact* fold angle (deviation from flatness)
    /// along each edge of the pattern, written as `(negative?, cos(a), sin(a))`.
    /// Both `cos(a)` and `sin(a)` are stored to ensure they're both in `N`.
    pub edges_fold_angle_cs: Option<TiVec<Edge, (bool, BasedExpr, BasedExpr)>>,

    /// For each edge, the approximate length of the edge.
    /// This is mainly useful for defining the intrinsic geometry of
    /// abstract complexes where `vertices_coords` are unspecified;
    /// otherwise, `edges_length` can be computed from `vertices_coords`.
    pub edges_length_f64: Option<TiVec<Edge, f64>>,
    /// For each edge, the exact squared length of the edge.
    /// This is mainly useful for defining the intrinsic geometry of
    /// abstract complexes where `vertices_coords` are unspecified;
    /// otherwise, `edges_length` can be computed from `vertices_coords`.
    pub edges_length2_exact: Option<TiVec<Edge, BasedExpr>>,

    /// For each face, an array *half-edge* IDs for the edges around
    /// the face *in counterclockwise order*. (See `HalfEdge`).
    /// 
    /// # Requirements
    /// * For each face, `edges_vertices[faces_half_edges[f][i]][1] == edges_vertices[faces_half_edges[f][(i+1)%n]][0]`.
    pub faces_half_edges: Option<TiVec<Face, Vec<HalfEdge>>>,
    
    /// An array of triples `(f, g, s)` where `f` and `g` are face IDs
    /// and `s` is a `FaceOrder`:
    /// * `Above` indicates that face `f` lies *above* face `g`,
    ///   i.e., on the side pointed to by `g`'s normal vector in the folded state.
    /// * `Below` indicates that face `f` lies *below* face `g`,
    ///   i.e., on the side opposite `g`'s normal vector in the folded state.
    /// * `Unknown` indicates that `f` and `g` have unknown stacking order
    ///   (e.g., they do not overlap in their interiors).
    ///
    /// **Recommended** for frames with interior-overlapping faces.
    pub face_orders: Option<Vec<(Face, Face, FaceOrder)>>,
    /// An array of triples `[e, f, s]` where `e` and `f` are edge IDs
    /// and `s` is a `EdgeOrder`.
    /// * `Left` indicates that edge `e` lies locally on the *left* side of edge `f`
    ///   (relative to edge `f`'s orientation given by `edges_vertices`)
    /// * `Right` indicates that edge `e` lies locally on the *right* side of edge
    ///   `f` (relative to edge `f`'s orientation given by `edges_vertices`)
    /// * 0 indicates that `e` and `f` have unknown stacking order
    ///   (e.g., they do not overlap in their interiors).
    ///
    /// This property makes sense only in 2D.
    /// **Recommended** for linkage configurations with interior-overlapping edges.
    pub edge_orders: Option<Vec<(Edge, Edge, EdgeOrder)>>,

    /// Parent frame ID.  Intuitively, this frame (the child)
    /// is a modification (or, in general, is related to) the parent frame.
    /// This property is optional, but enables organizing frames into a tree
    /// structure.
    pub frame_parent: Option<FrameIndex>,
    /// If true, any properties in the parent frame
    /// (or recursively inherited from an ancestor) that is not overridden in
    /// this frame are automatically inherited, allowing you to avoid duplicated
    /// data in many cases.
    pub frame_inherit: Option<bool>,
    /// Vertex custom data (has key `vertices_*` in file)
    pub vertices_custom: IndexMap<String, TiVec<Vertex, Value>>,
    /// Edge custom data (has key `edges_*` in file)
    pub edges_custom: IndexMap<String, TiVec<Edge, Value>>,
    /// Face custom data (has key `faces_*` in file)
    pub faces_custom: IndexMap<String, TiVec<Face, Value>>,
    /// Custom data
    pub other_custom: IndexMap<String, Value>,
}

impl Default for Frame {
    fn default() -> Self {
        Self {
            frame_author: Default::default(),
            frame_title: Default::default(),
            frame_description: Default::default(),
            frame_classes: Default::default(),
            frame_attributes: Default::default(),
            frame_unit: FrameUnit::Unit,
            basis: Default::default(),
            num_vertices: Default::default(),
            vertices_coords_f64: Default::default(),
            vertices_coords_exact: Default::default(),
            vertices_half_edges: Default::default(),
            edges_vertices: Default::default(),
            edges_faces: Default::default(),
            edges_assignment: Default::default(),
            edges_fold_angle_f64: Default::default(),
            edges_fold_angle_cs: Default::default(),
            edges_length_f64: Default::default(),
            edges_length2_exact: Default::default(),
            faces_half_edges: Default::default(),
            face_orders: Default::default(),
            edge_orders: Default::default(),
            frame_parent: Default::default(),
            frame_inherit: Default::default(),
            vertices_custom: Default::default(),
            edges_custom: Default::default(),
            faces_custom: Default::default(),
            other_custom: Default::default(),
        }
    }
}

#[cfg(test)]
mod test {
    use crate::fold::{EdgeAssignment, EdgeOrder, FaceOrder, FileClass, FrameAttribute, FrameClass, FrameUnit};

    #[test]
    fn test_serialize_file_classes() {
        let expecteds = ["\"singleModel\"", "\"multiModel\"", "\"animation\"", "\"diagrams\"", "\"test:crossSections\""];
        let inputs = [FileClass::SingleModel, FileClass::MultiModel, FileClass::Animation,
            FileClass::Diagrams, FileClass::Custom("test:crossSections".to_owned())];
        for (input, expected) in inputs.into_iter().zip(expecteds) {
            let output = serde_json::to_string(&input).unwrap_or_else(|e| panic!("{input:?} error: {e}"));
            assert_eq!(&output, expected);
        }
    }

    #[test]
    fn test_deserialize_file_classes() {
        let inputs = ["\"singleModel\"", "\"multiModel\"", "\"animation\"", "\"diagrams\"", "\"test:crossSections\""];
        let expecteds = [FileClass::SingleModel, FileClass::MultiModel, FileClass::Animation,
            FileClass::Diagrams, FileClass::Custom("test:crossSections".to_owned())];
        for (input, expected) in inputs.into_iter().zip(expecteds) {
            let output = serde_json::from_str::<FileClass>(input).unwrap_or_else(|e| panic!("{input} error: {e}"));
            assert_eq!(output, expected);
        }
    }

    #[test]
    fn test_serialize_frame_classes() {
        let expecteds = ["\"creasePattern\"", "\"foldedForm\"", "\"graph\"", "\"linkage\"", "\"test:triangleMesh\""];
        let inputs = [FrameClass::CreasePattern, FrameClass::FoldedForm, FrameClass::Graph,
            FrameClass::Linkage, FrameClass::Custom("test:triangleMesh".to_owned())];
        for (input, expected) in inputs.into_iter().zip(expecteds) {
            let output = serde_json::to_string(&input).unwrap_or_else(|e| panic!("{input:?} error: {e}"));
            assert_eq!(output, expected);
        }
    }

    #[test]
    fn test_deserialize_frame_classes() {
        let inputs = ["\"creasePattern\"", "\"foldedForm\"", "\"graph\"", "\"linkage\"", "\"test:triangleMesh\""];
        let expecteds = [FrameClass::CreasePattern, FrameClass::FoldedForm, FrameClass::Graph,
            FrameClass::Linkage, FrameClass::Custom("test:triangleMesh".to_owned())];
        for (input, expected) in inputs.into_iter().zip(expecteds) {
            let output = serde_json::from_str::<FrameClass>(input).unwrap_or_else(|e| panic!("{input} error: {e}"));
            assert_eq!(output, expected);
        }
    }

    #[test]
    fn test_serialize_frame_attributes() {
        let expecteds = ["\"2D\"", "\"3D\"", "\"abstract\"",
            "\"manifold\"", "\"nonManifold\"",
            "\"orientable\"", "\"nonOrientable\"",
            "\"selfTouching\"", "\"nonSelfTouching\"",
            "\"selfIntersecting\"", "\"nonSelfIntersecting\"",
            "\"cuts\"", "\"noCuts\"",
            "\"joins\"", "\"noJoins\"",
            "\"convexFaces\"", "\"nonConvexFaces\"",
            "\"test:torus\""];
        let inputs = [FrameAttribute::_2D, FrameAttribute::_3D, FrameAttribute::Abstract,
            FrameAttribute::Manifold, FrameAttribute::NonManifold,
            FrameAttribute::Orientable, FrameAttribute::NonOrientable,
            FrameAttribute::SelfTouching, FrameAttribute::NonSelfTouching,
            FrameAttribute::SelfIntersecting, FrameAttribute::NonSelfIntersecting,
            FrameAttribute::Cuts, FrameAttribute::NoCuts,
            FrameAttribute::Joins, FrameAttribute::NoJoins,
            FrameAttribute::ConvexFaces, FrameAttribute::NonConvexFaces,
            FrameAttribute::Custom("test:torus".to_owned())];
        for (input, expected) in inputs.into_iter().zip(expecteds) {
            let output = serde_json::to_string(&input).unwrap_or_else(|e| panic!("{input:?} error: {e}"));
            assert_eq!(output, expected);
        }
    }

    #[test]
    fn test_deserialize_frame_attributes() {
        let inputs = ["\"2D\"", "\"3D\"", "\"abstract\"",
            "\"manifold\"", "\"nonManifold\"",
            "\"orientable\"", "\"nonOrientable\"",
            "\"selfTouching\"", "\"nonSelfTouching\"",
            "\"selfIntersecting\"", "\"nonSelfIntersecting\"",
            "\"cuts\"", "\"noCuts\"",
            "\"joins\"", "\"noJoins\"",
            "\"convexFaces\"", "\"nonConvexFaces\"",
            "\"test:torus\""];
        let expecteds = [FrameAttribute::_2D, FrameAttribute::_3D, FrameAttribute::Abstract,
            FrameAttribute::Manifold, FrameAttribute::NonManifold,
            FrameAttribute::Orientable, FrameAttribute::NonOrientable,
            FrameAttribute::SelfTouching, FrameAttribute::NonSelfTouching,
            FrameAttribute::SelfIntersecting, FrameAttribute::NonSelfIntersecting,
            FrameAttribute::Cuts, FrameAttribute::NoCuts,
            FrameAttribute::Joins, FrameAttribute::NoJoins,
            FrameAttribute::ConvexFaces, FrameAttribute::NonConvexFaces,
            FrameAttribute::Custom("test:torus".to_owned())];
        for (input, expected) in inputs.into_iter().zip(expecteds) {
            let output = serde_json::from_str::<FrameAttribute>(input).unwrap_or_else(|e| panic!("{input} error: {e}"));
            assert_eq!(output, expected);
        }
    }

    #[test]
    fn test_serialize_frame_units() {
        let expecteds = ["\"unit\"", "\"in\"", "\"pt\"", "\"m\"", "\"cm\"", "\"mm\"", "\"um\"", "\"nm\"", "\"test:light-year\""];
        let inputs = [FrameUnit::Unit, FrameUnit::Inch, FrameUnit::Point,
            FrameUnit::Meter, FrameUnit::Centimeter, FrameUnit::Millimeter, FrameUnit::Micrometer, FrameUnit::Nanometer,
            FrameUnit::Custom("test:light-year".to_owned())];
        for (input, expected) in inputs.into_iter().zip(expecteds) {
            let output = serde_json::to_string(&input).unwrap_or_else(|e| panic!("{input:?} error: {e}"));
            assert_eq!(output, expected);
        }
    }

    #[test]
    fn test_deserialize_frame_units() {
        let inputs = ["\"unit\"", "\"in\"", "\"pt\"", "\"m\"", "\"cm\"", "\"mm\"", "\"um\"", "\"nm\"", "\"test:light-year\""];
        let expecteds = [FrameUnit::Unit, FrameUnit::Inch, FrameUnit::Point,
            FrameUnit::Meter, FrameUnit::Centimeter, FrameUnit::Millimeter, FrameUnit::Micrometer, FrameUnit::Nanometer,
            FrameUnit::Custom("test:light-year".to_owned())];
        for (input, expected) in inputs.into_iter().zip(expecteds) {
            let output = serde_json::from_str::<FrameUnit>(input).unwrap_or_else(|e| panic!("{input} error: {e}"));
            assert_eq!(output, expected);
        }
    }

    #[test]
    fn test_serialize_edge_assignments() {
        let expecteds = ["\"B\"", "\"F\"", "\"V\"", "\"M\"", "\"U\"", "\"C\"", "\"J\""];
        let inputs = [EdgeAssignment::Boundary, EdgeAssignment::Flat,
            EdgeAssignment::Valley, EdgeAssignment::Mountain,
            EdgeAssignment::Unassigned, EdgeAssignment::Cut, EdgeAssignment::Join];
        for (input, expected) in inputs.into_iter().zip(expecteds) {
            let output = serde_json::to_string(&input).unwrap_or_else(|e| panic!("{input:?} error: {e}"));
            assert_eq!(output, expected);
        }
    }

    #[test]
    fn test_deserialize_edge_assignments() {
        let inputs = ["\"B\"", "\"F\"", "\"V\"", "\"M\"", "\"U\"", "\"C\"", "\"J\""];
        let expecteds = [EdgeAssignment::Boundary, EdgeAssignment::Flat,
            EdgeAssignment::Valley, EdgeAssignment::Mountain,
            EdgeAssignment::Unassigned, EdgeAssignment::Cut, EdgeAssignment::Join];
        for (input, expected) in inputs.into_iter().zip(expecteds) {
            let output = serde_json::from_str::<EdgeAssignment>(input).unwrap_or_else(|e| panic!("{input} error: {e}"));
            assert_eq!(output, expected);
        }
    }

    #[test]
    fn test_serialize_face_orders() {
        let expecteds = ["-1", "0", "1"];
        let inputs = [FaceOrder::Below, FaceOrder::Unknown, FaceOrder::Above];
        for (input, expected) in inputs.into_iter().zip(expecteds) {
            let output = serde_json::to_string(&input).unwrap_or_else(|e| panic!("{input:?} error: {e}"));
            assert_eq!(output, expected);
        }
    }

    #[test]
    fn test_deserialize_face_orders() {
        let inputs = ["-1", "0", "1"];
        let expecteds = [FaceOrder::Below, FaceOrder::Unknown, FaceOrder::Above];
        for (input, expected) in inputs.into_iter().zip(expecteds) {
            let output = serde_json::from_str::<FaceOrder>(input).unwrap_or_else(|e| panic!("{input} error: {e}"));
            assert_eq!(output, expected);
        }
    }

    #[test]
    fn test_serialize_edge_orders() {
        let expecteds = ["-1", "0", "1"];
        let inputs = [EdgeOrder::Right, EdgeOrder::Unknown, EdgeOrder::Left];
        for (input, expected) in inputs.into_iter().zip(expecteds) {
            let output = serde_json::to_string(&input).unwrap_or_else(|e| panic!("{input:?} error: {e}"));
            assert_eq!(output, expected);
        }
    }

    #[test]
    fn test_deserialize_edge_orders() {
        let inputs = ["-1", "0", "1"];
        let expecteds = [EdgeOrder::Right, EdgeOrder::Unknown, EdgeOrder::Left];
        for (input, expected) in inputs.into_iter().zip(expecteds) {
            let output = serde_json::from_str::<EdgeOrder>(input).unwrap_or_else(|e| panic!("{input} error: {e}"));
            assert_eq!(output, expected);
        }
    }
}