use std::fmt;

// mostly copying the dbg! definition
#[macro_export]
macro_rules! my_dbg {
    () => {
        eprintln!("[{}:{}]", file!(), line!());
    };
    (inline $val:expr $(,)?) => {
        match $val {
            tmp => {
                eprintln!("[{}:{}] {} = {:?}",
                    file!(), line!(), stringify!($val), $crate::specialized_dbg::DebugWrapper(&tmp));
                tmp
            }
        }
    };
    ($val:expr $(,)?) => {
        match $val {
            tmp => {
                eprintln!("[{}:{}] {} = {:#?}",
                    file!(), line!(), stringify!($val), $crate::specialized_dbg::DebugWrapper(&tmp));
                tmp
            }
        }
    };
    (inline $($val:expr),+ $(,)?) => {
        ($(my_dbg!($val)),+,)
    };
    ($($val:expr),+ $(,)?) => {
        ($(my_dbg!(inline $val)),+,)
    };
}

pub(crate) struct DebugWrapper<'a, T: ?Sized>(pub(crate) &'a T);
impl<T: ?Sized> fmt::Debug for DebugWrapper<'_, T> {
    default fn fmt(&self, f: &mut fmt::Formatter<'_>) -> Result<(), fmt::Error> {
        write!(f, "{{{}}}", std::any::type_name::<T>())
    }
}
impl<T: ?Sized + fmt::Debug> fmt::Debug for DebugWrapper<'_, T> {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> Result<(), fmt::Error> {
        fmt::Debug::fmt(self.0, f)
    }
}