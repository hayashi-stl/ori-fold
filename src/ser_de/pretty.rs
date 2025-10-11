//! Pretty formatter for the FOLD format.
//! Differs from `serde_json`'s pretty formatter in that
//! elements of arrays are inlined.
//! Except `file_frames`.
//! 
//! So instead of 
//! ```json
//! "vertices_coords": [
//!     [
//!         0,
//!         1
//!     ],
//!     [
//!         2,
//!         3
//!     ]
//! ]
//! ```
//! we have
//! ```json
//! "vertices_coords": [
//!     [0, 1],
//!     [2, 3]
//! ]
//! ```
//! Much easier to read.
use std::io;

use serde::Serialize;
use serde_json::{ser::Formatter, Serializer};
use serde_json::Result;

/// Custom Result macro from `serde_json`
macro_rules! tri {
    ($e:expr $(,)?) => {
        match $e {
            core::result::Result::Ok(val) => val,
            core::result::Result::Err(err) => return core::result::Result::Err(err),
        }
    };
}

#[derive(Clone, Copy, Debug, PartialEq, Eq)]
enum FileFramesState {
    Outside,
    InObjectKey,
    Matched,
    Inside,
}

/// This structure pretty prints a JSON value to make it human readable.
#[derive(Clone, Debug)]
pub struct FoldFormatter<'a> {
    current_indent: usize,
    has_value: bool,
    array_level: usize,
    file_frames_state: FileFramesState,
    indent: &'a [u8],
}

impl<'a> FoldFormatter<'a> {
    /// Construct a pretty printer formatter that defaults to using four spaces for indentation.
    pub fn new() -> Self {
        FoldFormatter::with_indent(b"    ")
    }

    /// Construct a pretty printer formatter that uses the `indent` string for indentation.
    pub fn with_indent(indent: &'a [u8]) -> Self {
        FoldFormatter {
            current_indent: 0,
            has_value: false,
            array_level: 0,
            file_frames_state: FileFramesState::Outside,
            indent,
        }
    }

    fn should_indent_array(&self) -> bool {
        match self.file_frames_state {
            FileFramesState::Outside => self.array_level < 2,
            FileFramesState::InObjectKey => self.array_level < 2,
            FileFramesState::Matched => self.array_level < 2,
            FileFramesState::Inside => self.array_level < 3,
        }
    }
}

impl<'a> Default for FoldFormatter<'a> {
    fn default() -> Self {
        FoldFormatter::new()
    }
}

impl<'a> Formatter for FoldFormatter<'a> {
    #[inline]
    fn write_string_fragment<W>(&mut self, writer: &mut W, fragment: &str) -> io::Result<()> where W: ?Sized + io::Write {
        use FileFramesState::*;

        // key is not exactly "file_frames"
        if self.file_frames_state == Matched {
            self.file_frames_state = Outside;
        } else if self.file_frames_state == InObjectKey {
            self.file_frames_state = if fragment == "file_frames" { Matched } else { Outside };
        }
        writer.write_all(fragment.as_bytes())
    }

    fn write_char_escape<W>(&mut self, writer: &mut W, char_escape: serde_json::ser::CharEscape) -> io::Result<()> where W: ?Sized + io::Write, {
        use serde_json::ser::CharEscape::*;

        // key is not exactly "file_frames"
        if self.file_frames_state == FileFramesState::InObjectKey || self.file_frames_state == FileFramesState::Matched {
            self.file_frames_state = FileFramesState::Outside;
        }

        let escape_char = match char_escape {
            Quote => b'"',
            ReverseSolidus => b'\\',
            Solidus => b'/',
            Backspace => b'b',
            FormFeed => b'f',
            LineFeed => b'n',
            CarriageReturn => b'r',
            Tab => b't',
            AsciiControl(_) => b'u',
        };

        match char_escape {
            AsciiControl(byte) => {
                static HEX_DIGITS: [u8; 16] = *b"0123456789abcdef";
                let bytes = &[
                    b'\\',
                    escape_char,
                    b'0',
                    b'0',
                    HEX_DIGITS[(byte >> 4) as usize],
                    HEX_DIGITS[(byte & 0xF) as usize],
                ];
                writer.write_all(bytes)
            }
            _ => writer.write_all(&[b'\\', escape_char]),
        }
    }

    #[inline]
    fn begin_array<W>(&mut self, writer: &mut W) -> io::Result<()>
    where
        W: ?Sized + io::Write,
    {
        self.array_level += 1;
        if self.should_indent_array() {
            self.current_indent += 1;
        }
        self.has_value = false;
        writer.write_all(b"[")
    }

    #[inline]
    fn end_array<W>(&mut self, writer: &mut W) -> io::Result<()>
    where
        W: ?Sized + io::Write,
    {
        if self.should_indent_array() {
            self.current_indent -= 1;

            if self.has_value {
                tri!(writer.write_all(b"\n"));
                tri!(indent(writer, self.current_indent, self.indent));
            }
        }

        self.array_level -= 1;
        writer.write_all(b"]")
    }

    #[inline]
    fn begin_array_value<W>(&mut self, writer: &mut W, first: bool) -> io::Result<()>
    where
        W: ?Sized + io::Write,
    {
        if self.should_indent_array() {
            tri!(writer.write_all(if first { b"\n" } else { b",\n" }));
            tri!(indent(writer, self.current_indent, self.indent));
        } else {
            tri!(writer.write_all(if first { b"" } else { b", " }));
        }
        Ok(())
    }

    #[inline]
    fn end_array_value<W>(&mut self, _writer: &mut W) -> io::Result<()>
    where
        W: ?Sized + io::Write,
    {
        self.has_value = true;
        Ok(())
    }

    #[inline]
    fn begin_object<W>(&mut self, writer: &mut W) -> io::Result<()>
    where
        W: ?Sized + io::Write,
    {
        self.current_indent += 1;
        self.has_value = false;
        writer.write_all(b"{")
    }

    #[inline]
    fn end_object<W>(&mut self, writer: &mut W) -> io::Result<()>
    where
        W: ?Sized + io::Write,
    {
        self.current_indent -= 1;

        if self.has_value {
            tri!(writer.write_all(b"\n"));
            tri!(indent(writer, self.current_indent, self.indent));
        }

        writer.write_all(b"}")
    }

    #[inline]
    fn begin_object_key<W>(&mut self, writer: &mut W, first: bool) -> io::Result<()>
    where
        W: ?Sized + io::Write,
    {
        if self.file_frames_state == FileFramesState::Outside {
            self.file_frames_state = FileFramesState::InObjectKey;
        }
        tri!(writer.write_all(if first { b"\n" } else { b",\n" }));
        indent(writer, self.current_indent, self.indent)
    }

    #[inline]
    fn begin_object_value<W>(&mut self, writer: &mut W) -> io::Result<()>
    where
        W: ?Sized + io::Write,
    {
        if self.file_frames_state == FileFramesState::Matched {
            self.file_frames_state = FileFramesState::Inside;
        }
        writer.write_all(b": ")
    }

    #[inline]
    fn end_object_value<W>(&mut self, _writer: &mut W) -> io::Result<()>
    where
        W: ?Sized + io::Write,
    {
        if self.array_level == 0 {
            self.file_frames_state = FileFramesState::Outside;
        }
        self.has_value = true;
        Ok(())
    }
}

fn indent<W>(wr: &mut W, n: usize, s: &[u8]) -> io::Result<()>
where
    W: ?Sized + io::Write,
{
    for _ in 0..n {
        tri!(wr.write_all(s));
    }

    Ok(())
}

/// Serialize the given data structure as FOLD-printed JSON into the I/O
/// stream.
///
/// Serialization guarantees it only feeds valid UTF-8 sequences to the writer.
///
/// # Errors
///
/// Serialization can fail if `T`'s implementation of `Serialize` decides to
/// fail, or if `T` contains a map with non-string keys.
#[inline]
#[cfg_attr(docsrs, doc(cfg(feature = "std")))]
pub fn to_writer_fold<W, T>(writer: W, value: &T) -> Result<()>
where
    W: io::Write,
    T: ?Sized + Serialize,
{
    let mut ser = Serializer::with_formatter(writer, FoldFormatter::new());
    value.serialize(&mut ser)
}

/// Serialize the given data structure as a FOLD-printed JSON byte vector.
///
/// # Errors
///
/// Serialization can fail if `T`'s implementation of `Serialize` decides to
/// fail, or if `T` contains a map with non-string keys.
#[inline]
pub fn to_vec_fold<T>(value: &T) -> Result<Vec<u8>>
where
    T: ?Sized + Serialize,
{
    let mut writer = Vec::with_capacity(128);
    tri!(to_writer_fold(&mut writer, value));
    Ok(writer)
}

/// Serialize the given data structure as a FOLD-printed String of JSON.
///
/// # Errors
///
/// Serialization can fail if `T`'s implementation of `Serialize` decides to
/// fail, or if `T` contains a map with non-string keys.
#[inline]
pub fn to_string_fold<T>(value: &T) -> Result<String>
where
    T: ?Sized + Serialize,
{
    let vec = tri!(to_vec_fold(value));
    let string = unsafe {
        // We do not emit invalid UTF-8.
        String::from_utf8_unchecked(vec)
    };
    Ok(string)
}

