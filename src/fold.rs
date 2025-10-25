use std::{fmt::{Debug, Display}, fs::File, io::BufReader, ops::Index, path::Path};
use derive_more::{Add, AddAssign, Div, DivAssign, From, Into, Mul, MulAssign, Rem, RemAssign, Sub, SubAssign};
use exact_number::{basis::ArcBasis, BasedExpr, Angle};
use getset::{CopyGetters, Getters};
use indexmap::{IndexMap, IndexSet};
use nalgebra::{DVector, DMatrix};
use serde::{Deserialize, Serialize};
use serde_json::{Map, Value};
use serde_repr::{Deserialize_repr, Serialize_repr};
use typed_index_collections::{TiSlice, TiVec};

/// A subjective interpretation about what the entire file represents.
#[derive(Clone, Debug, PartialEq, Eq, PartialOrd, Ord, Hash)]
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
#[derive(Getters)]
#[cfg_attr(test, derive(PartialEq))] // Really no point in this derivation outside of tests
pub struct Fold {
    /// The version of the FOLD spec that the file assumes.
    /// **Strongly recommended**, in case we ever have to make
    /// backward-incompatible changes.
    #[getset(get = "pub")]
    pub(crate) file_spec: String,
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
    pub file_classes: IndexSet<FileClass>,
    /// The frames in the file. The key frame is frame 0.
    /// Do *not* empty the list; there must always be a frame 0.
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
#[derive(Clone, Debug, PartialEq, Eq, PartialOrd, Ord, Hash)]
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
#[derive(Clone, Debug, PartialEq, Eq, PartialOrd, Ord, Hash)]
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
    /// (this library treats this flag as just a note and does not check it)
    Abstract,
    /// the polyhedral complex is a manifold (has at most two faces incident to each edge)
    Manifold,
    /// the polyhedral complex is not a manifold
    /// (this library treats this flag as just a note and does not check it)
    NonManifold,
    /// the polyhedral complex is orientable, meaning it can be assigned a consistent normal direction (and hence it is also manifold)
    Orientable,
    /// the polyhedral complex is not orientable, meaning it cannot be assigned a consistent normal direction
    /// (this library treats this flag as just a note and does not check it)
    NonOrientable,
    /// the polyhedral complex has faces that touch in their relative interiors, in which case you probably want a face ordering
    /// (this library treats this flag as just a note and does not check it)
    SelfTouching,
    /// the polyhedral complex does not have faces that touch in their relative interiors
    NonSelfTouching,
    /// the polyhedral complex has properly intersecting faces
    /// (this library treats this flag as just a note and does not check it)
    SelfIntersecting,
    /// the polyhedral complex does not have properly intersecting faces
    NonSelfIntersecting,
    /// an edge has an assignment of [`Cut`](crate::EdgeAssignment::Cut) (cut/slit representing multiple [`Boundary`](crate::EdgeAssignment::Boundary) edges)
    /// (this library treats this flag as just a note and does not check it)
    Cuts,
    /// no edges have an assignment of [`Cut`](crate::EdgeAssignment::Cut) (cut/slit representing multiple [`Boundary`](crate::EdgeAssignment::Boundary) edges)
    NoCuts,
    /// an edge has an assignment of [`Join`](crate::EdgeAssignment::Boundary) (join)
    /// (this library treats this flag as just a note and does not check it)
    Joins,
    /// no edges have an assignment of [`Join`](crate::EdgeAssignment::Boundary) (join)
    NoJoins,
    /// all faces are convex polygons
    ConvexFaces,
    /// some faces are nonconvex
    /// (this library treats this flag as just a note and does not check it)
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
    /// cut/slit edge (should be treated as multiple [`Boundary`](crate::EdgeAssignment::Boundary) edges)
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

impl Edge {
    /// Split this edge into its half-edges
    pub fn split(self) -> [HalfEdge; 2] {
        [HalfEdge::new(self, false), HalfEdge::new(self, true)]
    }
}

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

/// Type-safe face corner index.
#[derive(Clone, Copy, PartialEq, Eq, Hash, PartialOrd, Ord)]
pub struct FaceCorner(pub Face, pub usize);

impl FaceCorner {
    pub fn new(face: Face, corner: usize) -> Self {
        Self(face, corner)
    }

    pub fn face(self) -> Face {
        self.0
    }

    pub fn corner(self) -> usize {
        self.1
    }

    /// Gets the previous corner in the face
    pub fn prev(self, faces: &FacesHalfEdgesSlice) -> Self {
        Self(self.0, (self.1 + faces[self.0].len() - 1) % faces[self.0].len())
    }

    /// Gets the next corner in the face
    pub fn next(self, faces: &FacesHalfEdgesSlice) -> Self {
        Self(self.0, (self.1 + 1) % faces[self.0].len())
    }
}

impl Debug for FaceCorner {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        write!(f, "f_{}c{}", self.0.0, self.1)
    }
}

impl Display for FaceCorner {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        <FaceCorner as Debug>::fmt(&self, f)
    }
}

/// Type-safe frame index.
#[derive(Clone, Copy, Debug, PartialEq, Eq, Hash, PartialOrd, Ord, From, Into,
    Add, AddAssign, Sub, SubAssign, Mul, MulAssign, Div, DivAssign, Rem, RemAssign)]
#[derive(Serialize, Deserialize)]
#[repr(transparent)]
#[serde(transparent)]
pub struct FrameIndex(pub usize);

pub type VerticesHalfEdges = TiVec<Vertex, Vec<HalfEdge>>;
pub type VerticesHalfEdgesSlice = TiSlice<Vertex, Vec<HalfEdge>>;
pub type EdgesVertices = TiVec<Edge, [Vertex; 2]>;
pub type EdgesVerticesSlice = TiSlice<Edge, [Vertex; 2]>;
pub type EdgesFaceCorners = TiVec<Edge, [Vec<FaceCorner>; 2]>;
pub type EdgesFaceCornersSlice = TiSlice<Edge, [Vec<FaceCorner>; 2]>;
pub type FacesHalfEdges = TiVec<Face, Vec<HalfEdge>>;
pub type FacesHalfEdgesSlice = TiSlice<Face, Vec<HalfEdge>>;
pub type VerticesCustom = TiVec<Vertex, Value>;
pub type VerticesCustomSlice = TiSlice<Vertex, Value>;
pub type EdgesCustom = TiVec<Edge, Value>;
pub type EdgesCustomSlice = TiSlice<Edge, Value>;
pub type FacesCustom = TiVec<Face, Value>;
pub type FacesCustomSlice = TiSlice<Face, Value>;

pub trait EdgesVerticesEx {
    fn at(&self, index: HalfEdge) -> [Vertex; 2];
    fn at_mut(&mut self, index: HalfEdge) -> [&mut Vertex; 2];
    fn try_at(&self, index: HalfEdge) -> Option<[Vertex; 2]>;
    fn try_at_mut(&mut self, index: HalfEdge) -> Option<[&mut Vertex; 2]>;
}

pub trait EdgesFaceCornersEx {
    fn at(&self, index: HalfEdge) -> &[FaceCorner];
    fn at_mut(&mut self, index: HalfEdge) -> &mut Vec<FaceCorner>;
    fn pair_at(&self, index: HalfEdge) -> [&[FaceCorner]; 2];
    fn try_at(&self, index: HalfEdge) -> Option<&[FaceCorner]>;
    fn try_at_mut(&mut self, index: HalfEdge) -> Option<&mut Vec<FaceCorner>>;
    fn try_pair_at(&self, index: HalfEdge) -> Option<[&[FaceCorner]; 2]>;
    fn half_iter_enumerated(&self) -> impl Iterator<Item = (HalfEdge, &Vec<FaceCorner>)>;
    fn half_iter_mut_enumerated(&mut self) -> impl Iterator<Item = (HalfEdge, &mut Vec<FaceCorner>)>;
}

impl EdgesVerticesEx for EdgesVerticesSlice {
    fn at(&self, index: HalfEdge) -> [Vertex; 2] {
        self.try_at(index).unwrap()
    }

    fn at_mut(&mut self, index: HalfEdge) -> [&mut Vertex; 2] {
        self.try_at_mut(index).unwrap()
    }

    fn try_at(&self, index: HalfEdge) -> Option<[Vertex; 2]> {
        let mut result = *self.get(index.edge())?;
        if index.flip_bit() { result.reverse() }
        Some(result)
    }

    fn try_at_mut(&mut self, index: HalfEdge) -> Option<[&mut Vertex; 2]> {
        let result = self.get_mut(index.edge())?;
        let mut result = result.get_disjoint_mut([0, 1]).unwrap();
        if index.flip_bit() { result.reverse() }
        Some(result)
    }
}

impl EdgesFaceCornersEx for EdgesFaceCornersSlice {
    fn at(&self, index: HalfEdge) -> &[FaceCorner] {
        self.try_at(index).unwrap()
    }

    fn at_mut(&mut self, index: HalfEdge) -> &mut Vec<FaceCorner> {
        self.try_at_mut(index).unwrap()
    }

    fn pair_at(&self, index: HalfEdge) -> [&[FaceCorner]; 2] {
        self.try_pair_at(index).unwrap()
    }

    fn try_at(&self, index: HalfEdge) -> Option<&[FaceCorner]> {
        self.get(index.edge())?.get(index.flip_bit() as usize).map(|cs| cs.as_slice())
    }

    fn try_at_mut(&mut self, index: HalfEdge) -> Option<&mut Vec<FaceCorner>> {
        self.get_mut(index.edge())?.get_mut(index.flip_bit() as usize)
    }

    fn try_pair_at(&self, index: HalfEdge) -> Option<[&[FaceCorner]; 2]> {
        self.try_at(index).and_then(|cs| self.try_at(index.flipped()).map(|cs2| [cs, cs2]))
    }

    fn half_iter_enumerated(&self) -> impl Iterator<Item = (HalfEdge, &Vec<FaceCorner>)> {
        self.iter_enumerated()
            .flat_map(|(e, cs)| [
                (HalfEdge::new(e, false), &cs[0]),
                (HalfEdge::new(e, true), &cs[1]),
            ])
    }

    fn half_iter_mut_enumerated(&mut self) -> impl Iterator<Item = (HalfEdge, &mut Vec<FaceCorner>)> {
        self.iter_mut_enumerated()
            .flat_map(|(e, cs)| {
                let (cs0, cs1) = cs.split_at_mut(1);
                [
                    (HalfEdge::new(e, false), &mut cs0[0]),
                    (HalfEdge::new(e, true), &mut cs1[0]),
                ]
            })
    }
}

pub trait AtFaceCorner {
    type Output;
    fn at(&self, index: FaceCorner) -> Self::Output;
    fn at_mut(&mut self, index: FaceCorner) -> &mut Self::Output;
    fn try_at(&self, index: FaceCorner) -> Option<Self::Output>;
    fn try_at_mut(&mut self, index: FaceCorner) -> Option<&mut Self::Output>;
}

impl AtFaceCorner for FacesHalfEdgesSlice {
    type Output = HalfEdge;

    fn at(&self, index: FaceCorner) -> Self::Output {
        self.try_at(index).unwrap()
    }

    fn at_mut(&mut self, index: FaceCorner) -> &mut Self::Output {
        self.try_at_mut(index).unwrap()
    }

    fn try_at(&self, index: FaceCorner) -> Option<Self::Output> {
        self.get(index.0)?.get(index.1).copied()
    }

    fn try_at_mut(&mut self, index: FaceCorner) -> Option<&mut Self::Output> {
        self.get_mut(index.0)?.get_mut(index.1)
    }
}

/// A type of vertex data
#[derive(Clone, Debug, PartialEq, Eq, Hash)]
pub enum VertexField {
    CoordsF64,
    CoordsExact,
    HalfEdges,
    Custom(String),
}

/// Vertex data. All data is documented in [`Frame`](crate::Frame)
#[derive(Clone, Debug, PartialEq)]
pub struct VertexData {
    pub coords_f64: Option<DVector<f64>>,
    pub coords_exact: Option<DVector<BasedExpr>>,
    pub half_edges: Option<Vec<HalfEdge>>,
    pub custom: IndexMap<String, Value>,
}

/// Vertex data for multiple vertices. All data is documented in [`Frame`](crate::Frame)
#[derive(Clone, Debug, PartialEq)]
pub struct VertexDatas {
    /// Needed in case all fields are undefined
    pub num_vertices: Option<usize>,
    pub coords_f64: Option<DMatrix<f64>>,
    pub coords_exact: Option<DMatrix<BasedExpr>>,
    pub half_edges: Option<Vec<Vec<HalfEdge>>>,
    pub custom: IndexMap<String, Vec<Value>>,
}

impl VertexDatas {
    /// Gets the number of elements whose data is stored.
    /// 
    /// # Panics
    /// Panics if there's disagreement.
    pub fn len(&self) -> usize {
        let counts = vec![
            self.num_vertices,
            self.coords_f64  .as_ref().map(|v| v.len()),
            self.coords_exact.as_ref().map(|v| v.len()),
            self.half_edges  .as_ref().map(|v| v.len()),
        ].into_iter().flatten().chain(
            self.custom.iter().map(|(_, v)| v.len())
        );
        counts.reduce(|curr, acc| if curr == acc { curr } else { panic!("disagreement over number of vertices in VertexDatas") })
            .expect("if no fields are defined, at least specify `num_vertices`")
    }
}

impl Default for VertexData {
    fn default() -> Self {
        Self {
            coords_f64: Default::default(),
            coords_exact: Default::default(),
            half_edges: Default::default(),
            custom: Default::default(),
        }
    }
}

impl Default for VertexDatas {
    fn default() -> Self {
        Self {
            num_vertices: Default::default(),
            coords_f64: Default::default(),
            coords_exact: Default::default(),
            half_edges: Default::default(),
            custom: Default::default(),
        }
    }
}

// A type of edge data
#[derive(Clone, Debug, PartialEq, Eq, Hash)]
pub enum EdgeField {
    Vertices,
    FaceCorners,
    Assignment,
    FoldAngleF64,
    FoldAngleExact,
    LengthF64,
    Length2Exact,
    Custom(String),
}


/// Edge data. All data is documented in [`Frame`](crate::Frame)
#[derive(Clone, Debug, PartialEq)]
pub struct EdgeData {
    pub vertices: [Vertex; 2],
    pub face_corners: Option<[Vec<FaceCorner>; 2]>,
    pub assignment: Option<EdgeAssignment>,
    pub fold_angle_f64: Option<f64>,
    pub fold_angle_exact: Option<Angle>,
    pub length_f64: Option<f64>,
    pub length2_exact: Option<BasedExpr>,
    pub custom: IndexMap<String, Value>,
}

/// Edge data for multiple edges. All data is documented in [`Frame`](crate::Frame)
#[derive(Clone, Debug, PartialEq)]
pub struct EdgeDatas {
    pub vertices: Vec<[Vertex; 2]>,
    pub face_corners: Option<Vec<[Vec<FaceCorner>; 2]>>,
    pub assignment: Option<Vec<EdgeAssignment>>,
    pub fold_angle_f64: Option<Vec<f64>>,
    pub fold_angle_exact: Option<Vec<Angle>>,
    pub length_f64: Option<Vec<f64>>,
    pub length2_exact: Option<Vec<BasedExpr>>,
    pub custom: IndexMap<String, Vec<Value>>,
}

impl EdgeData {
    /// `vertices` has no meaningful default,
    /// so just deal with this method
    pub fn default_with_vertices(vertices: [Vertex; 2]) -> Self {
        Self {
            vertices,
            face_corners: Default::default(),
            assignment: Default::default(),
            fold_angle_f64: Default::default(),
            fold_angle_exact: Default::default(),
            length_f64: Default::default(),
            length2_exact: Default::default(),
            custom: Default::default(),
        }
    }
}

impl EdgeDatas {
    /// `vertices` has no meaningful default,
    /// so just deal with this method
    pub fn default_with_vertices(vertices: Vec<[Vertex; 2]>) -> Self {
        Self {
            vertices,
            face_corners: Default::default(),
            assignment: Default::default(),
            fold_angle_f64: Default::default(),
            fold_angle_exact: Default::default(),
            length_f64: Default::default(),
            length2_exact: Default::default(),
            custom: Default::default(),
        }
    }
    
    /// Gets the number of elements whose data is stored.
    /// 
    /// # Panics
    /// Panics if there's disagreement.
    pub fn len(&self) -> usize {
        let counts = vec![
            Some(self.vertices.len()),
            self.face_corners    .as_ref().map(|v| v.len()),
            self.assignment      .as_ref().map(|v| v.len()),
            self.fold_angle_f64  .as_ref().map(|v| v.len()),
            self.fold_angle_exact.as_ref().map(|v| v.len()),
            self.length_f64      .as_ref().map(|v| v.len()),
            self.length2_exact   .as_ref().map(|v| v.len()),
        ].into_iter().flatten().chain(
            self.custom.iter().map(|(_, v)| v.len())
        );
        counts.reduce(|curr, acc| if curr == acc { curr } else { panic!("disagreement over number of edges in EdgeDatas") }).unwrap()
    }
}

/// A type of face data
#[derive(Clone, Debug, PartialEq, Eq, Hash)]
pub enum FaceField {
    HalfEdges,
    Custom(String),
}

/// Face data. All data is documented in [`Frame`](crate::Frame)
#[derive(Clone, Debug, PartialEq)]
pub struct FaceData {
    pub half_edges: Vec<HalfEdge>,
    pub custom: IndexMap<String, Value>,
}

/// Face data for multiple faces. All data is documented in [`Frame`](crate::Frame)
#[derive(Clone, Debug, PartialEq)]
pub struct FaceDatas {
    pub half_edges: Vec<Vec<HalfEdge>>,
    pub custom: IndexMap<String, Vec<Value>>,
}

impl FaceData {
    /// `half_edges` has no meaningful default,
    /// so just deal with this method
    pub fn default_with_half_edges(half_edges: Vec<HalfEdge>) -> Self {
        Self {
            half_edges,
            custom: Default::default(),
        }
    }
}

impl FaceDatas {
    /// `half_edges` has no meaningful default,
    /// so just deal with this method
    pub fn default_with_half_edges(half_edges: Vec<Vec<HalfEdge>>) -> Self {
        Self {
            half_edges,
            custom: Default::default(),
        }
    }

    /// Gets the number of elements whose data is stored.
    /// 
    /// # Panics
    /// Panics if there's disagreement.
    pub fn len(&self) -> usize {
        let counts = vec![
            Some(self.half_edges.len()),
        ].into_iter().flatten().chain(
            self.custom.iter().map(|(_, v)| v.len())
        );
        counts.reduce(|curr, acc| if curr == acc { curr } else { panic!("disagreement over number of faces in FaceDatas") }).unwrap()
    }
}

/// A FOLD frame, containing geometry information.
/// 
/// All operations prefer exact coordinates (`vertices_coords_exact`, `edges_fold_angle_exact`, `edges_length2_exact`),
/// not messing with approximate coordinates (`vertices_coords_f64`, `edges_fold_angle_f64`, `edges_length_f64`) when both exist.
#[derive(Clone, Debug)]
#[derive(Serialize, Deserialize)]
#[derive(Getters, CopyGetters)]
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
    pub frame_classes: IndexSet<FrameClass>,
    /// Attributes that objectively describe properties of the
    /// folded structure being represented.
    /// 
    /// # Requirements
    /// * At most 1 of [`_2D`](FrameAttribute::_2D), [`_3D`](FrameAttribute::_3D), and [`Abstract`](FrameAttribute::Abstract) are set.
    /// * [`Manifold`](FrameAttribute::Manifold) and [`NonManifold`](FrameAttribute::NonManifold) are not both set.
    /// * [`Orientable`](FrameAttribute::Orientable) and [`NonManifold`](FrameAttribute::NonManifold) are not both set.
    /// * [`Orientable`](FrameAttribute::Orientable) and [`NonOrientable`](FrameAttribute::NonOrientable) are not both set.
    /// * [`SelfTouching`](FrameAttribute::SelfTouching) and [`NonSelfTouching`](FrameAttribute::NonSelfTouching) are not both set.
    /// * [`SelfIntersecting`](FrameAttribute::SelfIntersecting) and [`NonSelfTouching`](FrameAttribute::NonSelfTouching) are not both set.
    /// * [`SelfIntersecting`](FrameAttribute::SelfIntersecting) and [`NonSelfIntersecting`](FrameAttribute::NonSelfIntersecting) are not both set.
    /// * [`Cuts`](FrameAttribute::Cuts) and [`NoCuts`](FrameAttribute::NoCuts) are not both set.
    /// * [`Joins`](FrameAttribute::Joins) and [`NoJoins`](FrameAttribute::NoJoins) are not both set.
    /// * [`ConvexFaces`](FrameAttribute::ConvexFaces) and [`NonConvexFaces`](FrameAttribute::NonConvexFaces) are not both set.
    #[getset(get = "pub")]
    pub(crate) frame_attributes: IndexSet<FrameAttribute>,
    /// Physical or logical unit that all coordinates are relative to
    pub frame_unit: FrameUnit,

    /// The basis that exact coordinates uses.
    #[getset(get = "pub")]
    pub(crate) basis: Option<ArcBasis>,

    /// For each vertex, an array of approximate coordinates,
    /// such as `[x, y, z]` or `[x, y]` (where `z` is implicitly zero).
    /// In higher dimensions, all trailing unspecified coordinates are implicitly
    /// zero.
    /// Vertex coordinates are columns.
    /// 
    /// Note that this isn't an array of `nalgebra::geometry::OPoint`
    /// because that one doesn't support a dynamic number of entries.
    /// 
    /// **Recommended** except for frames with attribute `Abstract`.
    /// 
    /// # Requirements
    /// * **Exists** if [`_2D`](FrameAttribute::_2D) or [`_3D`](FrameAttribute::_3D) is set
    ///   and `vertices_coords_exact == None`.
    /// * If [`_2D`](FrameAttribute::_2D) or [`_3D`](FrameAttribute::_3D) is set,
    ///   the coordinate dimensions match the attribute.
    #[getset(get = "pub")]
    pub(crate) vertices_coords_f64: Option<DMatrix<f64>>,
    /// For each vertex, an array of exact coordinates,
    /// such as `[x, y, z]` or `[x, y]` (where `z` is implicitly zero).
    /// In higher dimensions, all trailing unspecified coordinates are implicitly
    /// zero.
    /// Vertex coordinates are columns.
    /// 
    /// Note that this isn't an array of `nalgebra::geometry::OPoint`
    /// because that one doesn't support a dynamic number of entries.
    /// 
    /// **Recommended** except for frames with attribute `Abstract`.
    /// 
    /// # Requirements
    /// * If [`_2D`](FrameAttribute::_2D) or [`_3D`](FrameAttribute::_3D) is set,
    ///   the coordinate dimensions match the attribute.
    #[getset(get = "pub")]
    pub(crate) vertices_coords_exact: Option<DMatrix<BasedExpr>>,
    /// The number of vertices. Necessary to handle isolated vertices.
    #[getset(get_copy = "pub")]
    pub(crate) num_vertices: usize,
    /// For each vertex, an array of edge IDs for the *half-edge*s
    /// incident to the vertex.  If the frame represents an orientable manifold,
    /// this list should be ordered counterclockwise around the vertex.
    /// If the frame is a nonorientable manifold, this list should be cyclically
    /// ordered around the vertex.
    /// Note that `vertices_half_edges[v][i]` should be an edge *from* vertex
    /// `v` *to* the edge's other vertex.
    #[getset(get = "pub")]
    pub(crate) vertices_half_edges: Option<VerticesHalfEdges>,
    /// `edges_vertices`: For each edge, an array `[u, v]` of two different vertex IDs for
    /// the two endpoints of the edge.  This effectively defines the *orientation*
    /// of the edge, from `u` to `v`.  (This orientation choice is arbitrary,
    /// but is used to define the ordering of `edges_faces`.)
    /// **Recommended** in frames having any `edges_...` property
    /// (e.g., to represent mountain-valley assignment).
    #[getset(get = "pub")]
    pub(crate) edges_vertices: Option<EdgesVertices>,
    /// For each edge, an array of (face ID, corner index)s for the faces incident
    /// to the *unflipped half-edge* indicated by `edges_vertices`,
    /// followed by an array of (face ID, corner index)s for the faces incident
    /// to the *flipped half-edge*.
    /// For manifolds, the arrays for each edge should contain at most 2 face corners *total*.
    /// for orientable manifolds, the arrays for each edge should contain at most 1 face corner *each*.
    #[getset(get = "pub")]
    pub(crate) edges_face_corners: Option<EdgesFaceCorners>,
    /// For each edge, a string representing its fold direction assignment
    /// 
    /// This is *not* updated automatically when fold angles are updated.
    #[getset(get = "pub")]
    pub(crate) edges_assignment: Option<TiVec<Edge, EdgeAssignment>>,
    /// For each edge, the fold angle (deviation from flatness)
    /// along each edge of the pattern. The fold angle is positive for
    /// valley folds, negative for mountain folds, and zero for flat, unassigned,
    /// and border folds.  Accordingly, it is recommended for 
    /// the sign of `edge_foldAngle` to match `edges_assignment` if both are specified.
    /// 
    /// Note: this is *not* automatically updated when vertex coordinates or edge assignment are updated.
    #[getset(get = "pub")]
    pub(crate) edges_fold_angle_f64: Option<TiVec<Edge, f64>>,
    /// For each edge, the *exact* fold angle (deviation from flatness)
    /// along each edge of the pattern.
    /// 
    /// Note: this is *not* automatically updated when vertex coordinates or edge assignment are updated.
    #[getset(get = "pub")]
    pub(crate) edges_fold_angle_exact: Option<TiVec<Edge, Angle>>,

    /// For each edge, the approximate length of the edge.
    /// This is mainly useful for defining the intrinsic geometry of
    /// abstract complexes where `vertices_coords` are unspecified;
    /// otherwise, `edges_length` can be computed from `vertices_coords`.
    /// 
    /// Note: this is *not* automatically updated when vertex coordinates are updated.
    #[getset(get = "pub")]
    pub(crate) edges_length_f64: Option<TiVec<Edge, f64>>,
    /// For each edge, the exact squared length of the edge.
    /// This is mainly useful for defining the intrinsic geometry of
    /// abstract complexes where `vertices_coords` are unspecified;
    /// otherwise, `edges_length` can be computed from `vertices_coords`.
    /// 
    /// Note: this is *not* automatically updated when vertex coordinates are updated.
    #[getset(get = "pub")]
    pub(crate) edges_length2_exact: Option<TiVec<Edge, BasedExpr>>,

    /// For each face, an array *half-edge* IDs for the edges around
    /// the face *in counterclockwise order*. (See `HalfEdge`).
    #[getset(get = "pub")]
    pub(crate) faces_half_edges: Option<FacesHalfEdges>,
    
    /// An array of triples `(f, g, s)` where `f` and `g` are face IDs
    /// and `s` is a [`FaceOrder`]:
    /// * [`Above`](FaceOrder::Above) indicates that face `f` lies *above* face `g`,
    ///   i.e., on the side pointed to by `g`'s normal vector in the folded state.
    /// * [`Below`](FaceOrder::Below) indicates that face `f` lies *below* face `g`,
    ///   i.e., on the side opposite `g`'s normal vector in the folded state.
    /// * [`Unknown`](FaceOrder::Unknown) indicates that `f` and `g` have unknown stacking order
    ///   (e.g., they do not overlap in their interiors).
    ///
    /// **Recommended** for frames with interior-overlapping faces.
    #[getset(get = "pub")]
    pub(crate) face_orders: Option<Vec<(Face, Face, FaceOrder)>>,
    /// An array of triples `[e, f, s]` where `e` and `f` are edge IDs
    /// and `s` is a [`EdgeOrder`].
    /// * [`Left`](EdgeOrder::Left) indicates that edge `e` lies locally on the *left* side of edge `f`
    ///   (relative to edge `f`'s orientation given by `edges_vertices`)
    /// * [`Right`](EdgeOrder::Right) indicates that edge `e` lies locally on the *right* side of edge
    ///   `f` (relative to edge `f`'s orientation given by `edges_vertices`)
    /// * [`Unknown`](EdgeOrder::Unknown) indicates that `e` and `f` have unknown stacking order
    ///   (e.g., they do not overlap in their interiors).
    ///
    /// This property makes sense only in 2D.
    /// **Recommended** for linkage configurations with interior-overlapping edges.
    #[getset(get = "pub")]
    pub(crate) edge_orders: Option<Vec<(Edge, Edge, EdgeOrder)>>,

    /// Parent frame ID.  Intuitively, this frame (the child)
    /// is a modification (or, in general, is related to) the parent frame.
    /// This property is optional, but enables organizing frames into a tree
    /// structure.
    #[getset(get = "pub")]
    pub(crate) frame_parent: Option<FrameIndex>,
    /// If true, any properties in the parent frame
    /// (or recursively inherited from an ancestor) that is not overridden in
    /// this frame are automatically inherited, allowing you to avoid duplicated
    /// data in many cases.
    pub frame_inherit: Option<bool>,
    /// Vertex custom data (has key `vertices_*` in file)
    #[getset(get = "pub")]
    pub(crate) vertices_custom: IndexMap<String, TiVec<Vertex, Value>>,
    /// Edge custom data (has key `edges_*` in file)
    #[getset(get = "pub")]
    pub(crate) edges_custom: IndexMap<String, TiVec<Edge, Value>>,
    /// Face custom data (has key `faces_*` in file)
    #[getset(get = "pub")]
    pub(crate) faces_custom: IndexMap<String, TiVec<Face, Value>>,
    /// Custom data
    pub other_custom: IndexMap<String, Value>,
}

impl Frame {
    /// Adds an attribute without checking that the frame stays valid.
    /// Returns `true` if the attribute was not already there.
    pub(crate) fn add_attribute_unchecked(&mut self, attr: FrameAttribute) -> bool {
        self.frame_attributes.insert(attr)
    }

    /// Removes an attribute.
    /// Returns `true` if the attribute was there to remove.
    pub fn remove_attribute(&mut self, attr: FrameAttribute) -> bool {
        self.frame_attributes.swap_remove(&attr)
    }

    /// Gets a custom vertex field mutably. You cannot change the size of the array. 
    pub fn vertices_custom_field_mut(&mut self, field: &str) -> Option<&mut TiSlice<Vertex, Value>> {
        self.vertices_custom.get_mut(field).map(|v| v.as_mut_slice())
    }

    /// Gets a custom edge field mutably. You cannot change the size of the array. 
    pub fn edges_custom_field_mut(&mut self, field: &str) -> Option<&mut TiSlice<Edge, Value>> {
        self.edges_custom.get_mut(field).map(|v| v.as_mut_slice())
    }

    /// Gets a custom face field mutably. You cannot change the size of the array. 
    pub fn faces_custom_field_mut(&mut self, field: &str) -> Option<&mut TiSlice<Face, Value>> {
        self.faces_custom.get_mut(field).map(|v| v.as_mut_slice())
    }
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
            edges_face_corners: Default::default(),
            edges_assignment: Default::default(),
            edges_fold_angle_f64: Default::default(),
            edges_fold_angle_exact: Default::default(),
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