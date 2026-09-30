# Wrapping guide

This document attempts to outline some of the basic principles of the NekPy
wrapper, which relies entirely on the excellent `pybind11` binding library. An
extensive documentation is therefore beyond the scope of this document, but we
highlight aspects that are important for the NekPy wrappers.

In general, note that when things go wrong with `pybind11`, it'll be indicated
either by an extensive compiler error, or a runtime error in the Python
interpreter when you try to use your wrapper. Judicious use of Google is
therefore recommended to track down these issues!

You may also find the following resources useful:

- [The `pybind11` tutorial](https://pybind11.readthedocs.io/en/stable/basics.html)
- [The `pybind11` documentation](https://pybind11.readthedocs.io/en/stable/index.html)

# Basic structure

The NekPy wrapper is designed to mimic the library structure of Nektar++, with
directories for the `LibUtilities`, `SpatialDomains` and `StdRegions`
libraries. This is a deliberate design decision, so that classes and definitions
in Nektar++ can be easily located inside NekPy.

There are also some other directories and files:

- `NekPyConfig.hpp` is a convenience header that all `.cpp` files should
  import. It sets appropriate namespaces and macros that can be used for
  convenience.
- `lib` is a skeleton Python directory, which will be installed by CMake into
  the `dist` directory and automatically import the compiled binary libraries.
- `example` contains some basic Python examples.
- `cmake` has some CMake configuration files.

To demonstrate how to wrap classes, we'll refer to a number of existing parts of
the code below.

# Defining a library

First consider `LibUtilities`. An abbreviated version of the base file,
`LibUtilities.cpp` has the following structure:

```c++
#include <LibUtilities/Python/NekPyConfig.hpp>

void export_Basis(py::module &);
void export_SessionReader(py::module &);

// Define the LibUtilities module within the pybind11 module `m`.
PYBIND11_MODULE(_LibUtilities, m)
{
    // Export classes.
    export_Basis(m);
    export_SessionReader(m);
}
```

The `PYBIND11_MODULE(name, m)` macro allows us to define a Python module inside
C++. Note that in this case, the leading underscore in the name
(i.e. `_LibUtilities`) is deliberate. To define the contents of the module, we
call a number of functions that are prefixed by `export_`, which will define one
or more Python classes that live in this module. These Python classes correspond
with our Nektar++ classes. We adopt this approach simply so that we can split up
the different classes into different files, because it is not possible to call
`PYBIND11_MODULE` more than once. These functions are defined in appropriately
named files, for example:

- `export_Basis()` lives in the file
  `library/LibUtilities/Python/Foundations/Basis.cpp`
- This corresponds to the Nektar++ file
  `library/LibUtilities/Foundations/Basis.cpp` and the classes defined therein.

# Basic class wrapping

As a very basic example of wrapping a class, let's consider a minimal wrapper
for `SessionReader`.

```c++
void export_SessionReader(py::module &m)
{
    py::class_<SessionReader, std::shared_ptr<SessionReader>>(m,
                                                              "SessionReader")

        .def_static("CreateInstance", SessionReader_CreateInstance)

        .def("GetSessionName", &SessionReader::GetSessionName,
             py::return_value_policy::copy)

        .def("InitSession", &SessionReader::InitSession,
             py::arg("filenames") = py::list())

        .def("Finalise", &SessionReader::Finalise)
        ;
}
```

## `py::class_<>`

This `pybind11` object defines a Python class in C++. It is templated, and in
this case we have the following template arguments:

- `SessionReader` is the class that will be wrapped
- `std::shared_ptr<SessionReader>` indicates that this object should be stored
  inside a shared (or smart) pointer, which we frequently use throughout the
  library, as can be seen by the frequent use of `SessionReaderSharedPtr`

We then have two arguments to pass through:
- `m` is the module we are going to define the class within;
- `"SessionReader"` is the name of the class in Python.

## Wrapping member functions

We then call the `.def` function on the `class_<>`, which allows us to define
member functions on our class. This is equivalent to `def`-ing a function in
Python. `def` has two required parameters, and one optional parameter:

- The function name as a string, e.g. `"GetSessionName"`
- A function pointer that defines the C++ function that will be called
- An optional return policy, which we need to define when the C++ function
  returns a reference.

`pybind11` is very smart and can convert many Python objects to their equivalent
C++ function arguments, and C++ return types of the function to their respective
Python object. Many times therefore, one only needs to define the `.def()` call.

However, there are some instances where we need to do some additional
conversion, mask some C++ complexity from the Python interface, or deal with
functions that return references. We describe ways to deal with this below.

### Thin wrappers

Instead of defining a function pointer to a member of the C++ class, we can
define a function pointer to a separate function that defines some extra
functionality. This is called a **thin wrapper**.

As an example, consider the `CreateInstance` function. In C++ we pass this
function the command line arguments in the usual `argc`, `argv` format. In
Python, command line arguments are defined as a list of strings inside
`sys.argv`. However, `Boost.Python` does not know how to convert this list to
`argc, argv`, so we need some additional code.

```c++
SessionReaderSharedPtr SessionReader_CreateInstance(py::list &ns)
{
    // ... some code here that converts a Python list to the standard
    // c/c++ (int argc, char **argv) format for command line arguments.
    // Then use this to construct a SessionReader and return it.
    SessionReaderSharedPtr sr = SessionReader::CreateInstance(argc, argv);
    return sr;
}
```

In Python, we can then simply call 

```python
session = SessionReader.CreateInstance(sys.argv)
```

### Dealing with references

When dealing with functions in C++ that return references, e.g.

```c++
const NekDouble &GetFactor()
```

we need to supply an additional argument to `.def()`, since Python immutable
types such as strings and integers cannot be passed by reference. For a full
list of options, consult the `pybind11` [guide on return value
policies](https://pybind11.readthedocs.io/en/stable/advanced/functions.html#return-value-policies). In
the case above, we do not intend to change this from the Python side, and so we
can use `copy` as highlighted above, which will create a copy of the const
reference and return this to Python.

### Dealing with `Array<OneD, >`

The `LibUtilities/SharedArray.cpp` file contains a number of functions that
allow for the automatic conversion of Nektar++ `Array<OneD, >` storage to and
from NumPy `ndarray` objects. This means that you can wrap functions that take
these as parameters and return arrays very easily. However bear in mind the
following caveats:

- NumPy arrays created from Array objects will share their memory, so that
  changing the C++ array changes the contents of the NumPy array. Likewise, C++
  arrays created from NumPy arrays will share memory. Reference counting and
  capsules are used to ensure that memory should be persistent whilst the arrays
  are in use, either within Python or C++.
- Many functions in Nektar++ return Arrays through argument parameters. In
  Python this is a very unnatural way to write functions. For example:
  ```python
  # This is good
  x, y, z = exp.GetCoords()
  # This is bad
  x, y, z = np.zeros(10), np.zeros(10), np.zeros(10)
  exp.GetCoords(x,y,z)
  # Here I would assume x, y and z are updated as a result of GetCoords().
  ```
  Use thin wrappers to overcome this problem. For examples of how to do this,
  particularly in returning tuples, consult the `StdRegions/StdExpansion.cpp`
  wrapper.
- `TwoD` and `ThreeD` arrays are not presently supported.

## Inheritance

Nektar++ makes heavy use of inheritance, which can be translated to Python quite
easily using `pybind11`. For a good example of how to do this, you can examine
the `StdRegions` wrapper for `StdExpansion` and its elements such as
`StdQuadExp`. In a cut-down form, these look like the following:

```c++
void export_StdExpansion(py::module &m)
{
    py::class_<StdExpansion, std::shared_ptr<StdExpansion>>(m, "StdExpansion");
}
void export_StdQuadExp(py::module &m)
{
    py::class_<StdQuadExp, StdExpansion, std::shared_ptr<StdQuadExp>>(
        m, "StdQuadExp", py::multiple_inheritance())
        .def(py::init<const LibUtilities::BasisKey &,
                      const LibUtilities::BasisKey &>());
}
```

Note the following:
- `StdExpansion` is an abstract class, so it has no initialiser defined within
  the body of the `.def()` calls.
- We use some, but not all of the subclasses in the definition of `StdQuadExp`
  to define its parent class. This does not necessarily need to include the full
  hierarchy of C++ inheritance: in `StdRegions` the inheritance graph for
  `StdQuadExp` looks like
  ```
  StdExpansion -> StdExpansion2D -> StdQuadExp
  ```
  In the above wrapper, we omit the StdExpansion2D call entirely.
- However, in order to do this, we need to use the `py::multiple_inheritance()`
  as an additional argument to `py::class_`.
- `py::init<>` is used to show how to wrap a C++ constructor. This can accept
  any arguments for which you have either written explicit wrappers or
  `pybind11` already knows how to convert.

## Wrapping enums

Most Nektar++ enumerators come in the form:

```c++
enum MyEnum {
    eItemOne,
    eItemTwo,
    SIZE_MyEnum
};
static const char *MyEnumMap[] = {
    "ItemOne"
    "ItemTwo"
    "ItemThree"
};
```

To wrap this, you can use the `NEKPY_WRAP_ENUM` macro defined in
`NekPyConfig.hpp`, which in this case can be used as `NEKPY_WRAP_ENUM(m, MyEnum,
MyEnumMap)`. Note that if instead of `const char *` the map is defined as a
`const std::string`, you can use the `NEKPY_WRAP_ENUM_STRING` macro.
