///////////////////////////////////////////////////////////////////////////////
//
// File: H5.h
//
// For more information, please see: http://www.nektar.info
//
// The MIT License
//
// Copyright (c) 2006 Division of Applied Mathematics, Brown University (USA),
// Department of Aeronautics, Imperial College London (UK), and Scientific
// Computing and Imaging Institute, University of Utah (USA).
//
// Permission is hereby granted, free of charge, to any person obtaining a
// copy of this software and associated documentation files (the "Software"),
// to deal in the Software without restriction, including without limitation
// the rights to use, copy, modify, merge, publish, distribute, sublicense,
// and/or sell copies of the Software, and to permit persons to whom the
// Software is furnished to do so, subject to the following conditions:
//
// The above copyright notice and this permission notice shall be included
// in all copies or substantial portions of the Software.
//
// THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS
// OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
// FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL
// THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
// LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING
// FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER
// DEALINGS IN THE SOFTWARE.
//
// Description: Simple OO wrapper around HDF5
//
///////////////////////////////////////////////////////////////////////////////

#ifndef NEKTAR_LIB_UTILITIES_BASIC_UTILS_H5_H
#define NEKTAR_LIB_UTILITIES_BASIC_UTILS_H5_H

#include <exception>
#include <hdf5.h>
#include <memory>
#include <string>
#include <vector>

#include <LibUtilities/BasicUtils/ErrorUtil.hpp>
#include <LibUtilities/Communication/Comm.h>

namespace Nektar::LibUtilities::H5
{

#define H5_CONSTRUCT(ans, func, args)                                          \
    {                                                                          \
        hid_t ret = func args;                                                 \
        if (ret < 0)                                                           \
            ErrorUtil::Error(ErrorUtil::efatal, __FILE__, __LINE__,            \
                             "HDF5 error in API function " #func, 0);          \
        ans = ret;                                                             \
    }
#define H5_CALL(func, args)                                                    \
    {                                                                          \
        herr_t ret = func args;                                                \
        if (ret < 0)                                                           \
            ErrorUtil::Error(ErrorUtil::efatal, __FILE__, __LINE__,            \
                             "HDF5 error in API function " #func, 0);          \
    }

class Error : public std::exception
{
};

// Forward declare
class Object;
typedef std::shared_ptr<Object> ObjectSharedPtr;
class DataType;
typedef std::shared_ptr<DataType> DataTypeSharedPtr;
class CompoundDataType;
typedef std::shared_ptr<CompoundDataType> CompoundDataTypeSharedPtr;
class DataSpace;
typedef std::shared_ptr<DataSpace> DataSpaceSharedPtr;
class CanHaveAttributes;
typedef std::shared_ptr<CanHaveAttributes> CanHaveAttributesSharedPtr;
class Attribute;
typedef std::shared_ptr<Attribute> AttributeSharedPtr;
class CanHaveGroupsDataSets;
typedef std::shared_ptr<CanHaveGroupsDataSets> CanHaveGroupsDataSetsSharedPtr;
class Group;
typedef std::shared_ptr<Group> GroupSharedPtr;
class File;
typedef std::shared_ptr<File> FileSharedPtr;
class DataSet;
typedef std::shared_ptr<DataSet> DataSetSharedPtr;
class PList;
typedef std::shared_ptr<PList> PListSharedPtr;

/// HDF5 base class
class Object : public std::enable_shared_from_this<Object>
{
public:
    void Close()
    {
        v_Close();
    }
    inline hid_t GetId() const
    {
        return m_Id;
    }
    // Overload cast to the HDF5 ID type to make objects
    // transparently usable in the HDF5 API.
    inline operator hid_t() const
    {
        return GetId();
    }

protected:
    LIB_UTILITIES_EXPORT Object();
    LIB_UTILITIES_EXPORT Object(hid_t id);
    LIB_UTILITIES_EXPORT virtual ~Object();
    hid_t m_Id;
    LIB_UTILITIES_EXPORT virtual void v_Close() = 0;
};

// PropertyList objects
class PList : public Object
{
public:
    /// Default options
    LIB_UTILITIES_EXPORT static PListSharedPtr Default();
    /// Properties for object creation
    LIB_UTILITIES_EXPORT static PListSharedPtr ObjectCreate();
    /// Properties for file creation
    LIB_UTILITIES_EXPORT static PListSharedPtr FileCreate();
    /// Properties for file access
    LIB_UTILITIES_EXPORT static PListSharedPtr FileAccess();
    /// Properties for dataset creation
    LIB_UTILITIES_EXPORT static PListSharedPtr DatasetCreate();
    /// Properties for dataset access
    LIB_UTILITIES_EXPORT static PListSharedPtr DatasetAccess();
    /// Properties for raw data transfer
    LIB_UTILITIES_EXPORT static PListSharedPtr DatasetXfer();
    /// Properties for file mounting
    LIB_UTILITIES_EXPORT static PListSharedPtr FileMount();
    /// Properties for group creation
    LIB_UTILITIES_EXPORT static PListSharedPtr GroupCreate();
    /// Properties for group access
    LIB_UTILITIES_EXPORT static PListSharedPtr GroupAccess();
    /// Properties for datatype creation
    LIB_UTILITIES_EXPORT static PListSharedPtr DatatypeCreate();
    /// Properties for datatype access
    LIB_UTILITIES_EXPORT static PListSharedPtr DatatypeAccess();
    /// Properties for character encoding when encoding strings or object names
    LIB_UTILITIES_EXPORT static PListSharedPtr StringCreate();
    /// Properties for attribute creation
    LIB_UTILITIES_EXPORT static PListSharedPtr AttributeCreate();
    /// Properties governing the object copying process
    LIB_UTILITIES_EXPORT static PListSharedPtr ObjectCopy();
    /// Properties governing link creation
    LIB_UTILITIES_EXPORT static PListSharedPtr LinkCreate();
    /// Properties governing link traversal when accessing objects
    LIB_UTILITIES_EXPORT static PListSharedPtr LinkAccess();

    LIB_UTILITIES_EXPORT PList();
    LIB_UTILITIES_EXPORT ~PList() override;
    LIB_UTILITIES_EXPORT void SetChunk(const std::vector<hsize_t> &dims);
    LIB_UTILITIES_EXPORT void SetDeflate(const unsigned level = 1);
    LIB_UTILITIES_EXPORT void SetMpio(CommSharedPtr comm);
    LIB_UTILITIES_EXPORT void SetDxMpioCollective();
    LIB_UTILITIES_EXPORT void SetDxMpioIndependent();

protected:
    LIB_UTILITIES_EXPORT void v_Close() override;

private:
    PList(hid_t cls);
};

/// Mixin for objects that contain groups and datasets (Group and File)
class CanHaveGroupsDataSets : public virtual Object
{
public:
    class LinkIterator
    {
    public:
        LIB_UTILITIES_EXPORT LinkIterator(CanHaveGroupsDataSetsSharedPtr grp,
                                          hsize_t idx = 0);
        LIB_UTILITIES_EXPORT const std::string &operator*();
        LIB_UTILITIES_EXPORT LinkIterator &operator++();
        LIB_UTILITIES_EXPORT bool operator==(const LinkIterator &other) const;
        inline bool operator!=(const LinkIterator &other) const
        {
            return !(*this == other);
        }
        inline hsize_t GetPos() const
        {
            return m_idx;
        }

        inline std::string GetName() const
        {
            return m_currentName;
        }

    private:
        LIB_UTILITIES_EXPORT static herr_t helper(hid_t g_id, const char *name,
                                                  const H5L_info_t *info,
                                                  void *op_data);
        CanHaveGroupsDataSetsSharedPtr m_grp;
        hsize_t m_idx;
        hsize_t m_next;
        hsize_t m_size;
        std::string m_currentName;
    };

    // Create a group with the given name. The createPL can be
    // omitted to use the default properties.
    LIB_UTILITIES_EXPORT GroupSharedPtr CreateGroup(
        const std::string &name, PListSharedPtr createPL = PList::Default(),
        PListSharedPtr accessPL = PList::Default());

    // Create a dataset with the name, type and space.
    // The createPL can be omitted to use the defaults.
    LIB_UTILITIES_EXPORT DataSetSharedPtr CreateDataSet(
        const std::string &name, DataTypeSharedPtr type,
        DataSpaceSharedPtr space, PListSharedPtr createPL = PList::Default(),
        PListSharedPtr accessPL = PList::Default());

    // Create a dataset containing the data supplied
    // The createPL can be omitted to use the defaults
    template <class T>
    DataSetSharedPtr CreateWriteDataSet(
        const std::string &name, const std::vector<T> &data,
        PListSharedPtr createPL = PList::Default(),
        PListSharedPtr accessPL = PList::Default());

    // Open an existing group.
    // The accessPL can be omitted to use the defaults
    LIB_UTILITIES_EXPORT GroupSharedPtr
    OpenGroup(const std::string &name,
              PListSharedPtr accessPL = PList::Default()) const;

    // Open an existing dataset
    // The accessPL can be omitted to use the defaults
    LIB_UTILITIES_EXPORT DataSetSharedPtr
    OpenDataSet(const std::string &name,
                PListSharedPtr accessPL = PList::Default()) const;

    LIB_UTILITIES_EXPORT bool ContainsDataSet(std::string nm);

    hsize_t GetNumElements()
    {
        return v_GetNumElements();
    }

protected:
    LIB_UTILITIES_EXPORT virtual hsize_t v_GetNumElements() = 0;
    LinkIterator begin();
    LinkIterator end();

    friend class key_iterator;
};

/// Mixin for objects that can have attributes (Group, DataSet, DataType)
class CanHaveAttributes : public virtual Object
{
public:
    class AttrIterator
    {
    public:
        LIB_UTILITIES_EXPORT AttrIterator(CanHaveAttributesSharedPtr obj,
                                          hsize_t idx = 0);
        LIB_UTILITIES_EXPORT const std::string &operator*();
        LIB_UTILITIES_EXPORT AttrIterator &operator++();
        LIB_UTILITIES_EXPORT bool operator==(const AttrIterator &other) const;
        inline bool operator!=(const AttrIterator &other) const
        {
            return !(*this == other);
        }
        inline hsize_t GetPos() const
        {
            return m_idx;
        }

    private:
        static herr_t helper(hid_t g_id, const char *name,
                             const H5A_info_t *info, void *op_data);
        CanHaveAttributesSharedPtr m_obj;
        hsize_t m_idx;
        hsize_t m_next;
        hsize_t m_size;
        std::string m_currentName;
    };

    LIB_UTILITIES_EXPORT AttributeSharedPtr
    CreateAttribute(const std::string &name, DataTypeSharedPtr type,
                    DataSpaceSharedPtr space);
    LIB_UTILITIES_EXPORT AttributeSharedPtr
    OpenAttribute(const std::string &name);

    template <class T>
    void SetAttribute(const std::string &name, const T &value);
    template <class T>
    void SetAttribute(const std::string &name, const std::vector<T> &value);

    template <class T> void GetAttribute(const std::string &name, T &value);
    template <class T>
    void GetAttribute(const std::string &name, std::vector<T> &value);

    LIB_UTILITIES_EXPORT int GetNumAttr() const;
    LIB_UTILITIES_EXPORT AttrIterator attr_begin();
    LIB_UTILITIES_EXPORT AttrIterator attr_end();
};

/// HDF5 DataSpace wrapper
class DataSpace : public Object
{
public:
    LIB_UTILITIES_EXPORT static DataSpaceSharedPtr Null();
    LIB_UTILITIES_EXPORT static DataSpaceSharedPtr Scalar();
    LIB_UTILITIES_EXPORT static DataSpaceSharedPtr OneD(hsize_t size);

    LIB_UTILITIES_EXPORT DataSpace();
    LIB_UTILITIES_EXPORT DataSpace(hsize_t size,
                                   hsize_t max = H5S_UNLIMITED - 1);
    LIB_UTILITIES_EXPORT DataSpace(const std::vector<hsize_t> &dims);
    LIB_UTILITIES_EXPORT DataSpace(const std::vector<hsize_t> &dims,
                                   const std::vector<hsize_t> &max_dims);
    LIB_UTILITIES_EXPORT ~DataSpace() override;

    LIB_UTILITIES_EXPORT void SelectRange(const hsize_t start,
                                          const hsize_t count);
    LIB_UTILITIES_EXPORT void AppendRange(const hsize_t start,
                                          const hsize_t count);

    LIB_UTILITIES_EXPORT void SelectRange(const std::vector<hsize_t> start,
                                          const std::vector<hsize_t> count);
    LIB_UTILITIES_EXPORT void AppendRange(const std::vector<hsize_t> start,
                                          const std::vector<hsize_t> count);

    LIB_UTILITIES_EXPORT void SetSelection(const hsize_t num_elmt,
                                           const std::vector<hsize_t> &coords);

    LIB_UTILITIES_EXPORT void ClearRange();

    LIB_UTILITIES_EXPORT hsize_t GetSize();
    LIB_UTILITIES_EXPORT std::vector<hsize_t> GetDims();

protected:
    LIB_UTILITIES_EXPORT void v_Close() override;

private:
    DataSpace(hid_t id);
    friend class Attribute;
    friend class DataSet;
};

// Policy class for the DataTypesTraits controlling whether data
// has to be converted in anyway before writing. (It could perhaps
// be DataTypeTraitsTraits but that's too horrible a name.)
//
// Default policy is to not convert at all.
template <class T> struct DataTypeConversionPolicy
{
    static const bool MustConvert = false;
    typedef T ConvertedType;
    typedef T ConvertedVectorElemType;
    static ConvertedType Convert(const T &obj);
    static T Deconvert(const ConvertedType &obj);
};

/// Traits class for HDF5 data types.
template <class T> struct DataTypeTraits
{
    typedef DataTypeConversionPolicy<T> Converter;

    /***
     * Define this for a specialision for any HDF5 NATIVE type you want to use.
     * See http://hdfgroup.org/HDF5/doc/UG/UG_frame11Datatypes.html
     */
    LIB_UTILITIES_EXPORT static const hid_t NativeType;

    typedef typename Converter::ConvertedType ConvertedType;

    static ConvertedType Convert(const T &obj);
    static T Deconvert(const ConvertedType &obj);
    /**
     * Get the address of the start of the data.
     * Default implementation just uses "&"
     */
    static const void *GetAddress(const ConvertedType &obj);
    static void *GetAddress(ConvertedType &obj);
    /**
     * Return a DataType object representing T.
     * Default implementation just calls PredefinedDataType::Native<T>()
     */
    static DataTypeSharedPtr GetType();
    static DataTypeSharedPtr GetType(const T &obj);
};

/// Wrap and HDF5 data type object. Technically this can have attributes, but
/// not really bothered.
class DataType : public Object
{
public:
    LIB_UTILITIES_EXPORT static DataTypeSharedPtr String(size_t len = 0);
    template <class T>
    static DataTypeSharedPtr OfObject([[maybe_unused]] const T &obj)
    {
        return DataTypeTraits<T>::GetType();
    }
    LIB_UTILITIES_EXPORT DataTypeSharedPtr Copy() const;

protected:
    LIB_UTILITIES_EXPORT void v_Close() override;
    LIB_UTILITIES_EXPORT DataType(hid_t id);
};

class CompoundDataType : public DataType
{
public:
    LIB_UTILITIES_EXPORT static CompoundDataTypeSharedPtr Create(size_t sz);

    void Add(std::string name, size_t offset, hid_t type)
    {
        H5Tinsert(m_Id, name.c_str(), offset, type);
    }
    void AddString(std::string name, size_t offset, size_t size)
    {
        hid_t strtype = H5Tcopy(H5T_C_S1);
        H5Tset_size(strtype, size);
        H5Tinsert(m_Id, name.c_str(), offset, strtype);
    }

protected:
    LIB_UTILITIES_EXPORT void v_Close() override;

private:
    CompoundDataType(hid_t);
};

/// Predefined HDF data types that must not be closed when done with.
class PredefinedDataType : public DataType
{
public:
    template <class T> static DataTypeSharedPtr Native();
    LIB_UTILITIES_EXPORT static DataTypeSharedPtr CS1();

protected:
    LIB_UTILITIES_EXPORT void v_Close() override;

private:
    LIB_UTILITIES_EXPORT PredefinedDataType(hid_t);
};

/// HDF5 Attribute Wrapper
class Attribute : public Object
{
public:
    LIB_UTILITIES_EXPORT ~Attribute() override;
    LIB_UTILITIES_EXPORT DataSpaceSharedPtr GetSpace() const;

protected:
    LIB_UTILITIES_EXPORT void v_Close() override;

private:
    Attribute(hid_t id) : Object(id)
    {
    }
    static AttributeSharedPtr Create(hid_t parent, const std::string &name,
                                     DataTypeSharedPtr type,
                                     DataSpaceSharedPtr space);
    static AttributeSharedPtr Open(hid_t parent, const std::string &name);
    friend class CanHaveAttributes;
};

/// HDF5 file wrapper
class File : public CanHaveGroupsDataSets
{
public:
    LIB_UTILITIES_EXPORT static FileSharedPtr Create(
        const std::string &filename, unsigned mode,
        PListSharedPtr createPL = PList::Default(),
        PListSharedPtr accessPL = PList::Default());
    LIB_UTILITIES_EXPORT static FileSharedPtr Open(
        const std::string &filename, unsigned mode,
        PListSharedPtr accessPL = PList::Default());
    LIB_UTILITIES_EXPORT ~File() override;

protected:
    LIB_UTILITIES_EXPORT void v_Close() override;
    LIB_UTILITIES_EXPORT hsize_t v_GetNumElements() override;

private:
    File(hid_t id);
};

/// HDF5 Group wrapper
class Group : public CanHaveAttributes, public CanHaveGroupsDataSets
{
public:
    LIB_UTILITIES_EXPORT ~Group() override;
    LIB_UTILITIES_EXPORT std::vector<std::string> GetElementNames();
    LIB_UTILITIES_EXPORT CanHaveAttributesSharedPtr operator[](hsize_t idx);
    LIB_UTILITIES_EXPORT CanHaveAttributesSharedPtr
    operator[](const std::string &key);

protected:
    LIB_UTILITIES_EXPORT void v_Close() override;
    LIB_UTILITIES_EXPORT hsize_t v_GetNumElements() override;

private:
    Group(hid_t id);
    friend class CanHaveGroupsDataSets;
};

class DataSet : public CanHaveAttributes
{
public:
    LIB_UTILITIES_EXPORT ~DataSet() override;
    LIB_UTILITIES_EXPORT DataSpaceSharedPtr GetSpace() const;

    template <class T> void Write(const std::vector<T> &data)
    {
        DataTypeSharedPtr mem_t = DataTypeTraits<T>::GetType();
        H5_CALL(H5Dwrite, (m_Id, mem_t->GetId(), H5S_ALL, H5S_ALL, H5P_DEFAULT,
                           &data[0]));
    }
    template <class T>
    void Write(const std::vector<T> &data, DataSpaceSharedPtr filespace,
               PListSharedPtr dxpl = PList::Default())
    {
        DataTypeSharedPtr mem_t     = DataTypeTraits<T>::GetType();
        DataSpaceSharedPtr memspace = DataSpace::OneD(data.size());

        H5_CALL(H5Dwrite, (m_Id, mem_t->GetId(), memspace->GetId(),
                           filespace->GetId(), dxpl->GetId(), &data[0]));
    }

    template <class T>
    void Write(const std::vector<T> &data, DataSpaceSharedPtr filespace,
               DataTypeSharedPtr dt, PListSharedPtr dxpl = PList::Default())
    {
        DataSpaceSharedPtr memspace = DataSpace::OneD(data.size());

        H5Dwrite(m_Id, dt->GetId(), memspace->GetId(), filespace->GetId(),
                 dxpl->GetId(), &data[0]);
    }

    void WriteString(std::string s, DataSpaceSharedPtr filespace,
                     DataTypeSharedPtr type,
                     PListSharedPtr dxpl = PList::Default())
    {
        H5_CALL(H5Dwrite, (m_Id, type->GetId(), H5S_ALL, filespace->GetId(),
                           dxpl->GetId(), s.c_str()));
    }

    void WriteVectorString(std::vector<std::string> s,
                           DataSpaceSharedPtr filespace, DataTypeSharedPtr type,
                           PListSharedPtr dxpl = PList::Default())
    {
        const char **ret;
        ret = new const char *[s.size()];
        for (size_t i = 0; i < s.size(); ++i)
        {
            ret[i] = s[i].c_str();
        }
        H5Dwrite(m_Id, type->GetId(), H5S_ALL, filespace->GetId(),
                 dxpl->GetId(), ret);
    }

    template <class T> void Read(std::vector<T> &data)
    {
        DataTypeSharedPtr mem_t  = DataTypeTraits<T>::GetType();
        DataSpaceSharedPtr space = GetSpace();
        ASSERTL0(H5Sget_simple_extent_ndims(space->GetId()) == 1,
                 "vector data not 1D");
        hsize_t len, maxdim;
        H5Sget_simple_extent_dims(space->GetId(), &len, &maxdim);

        data.resize(len);

        H5_CALL(H5Dread, (m_Id, mem_t->GetId(), H5S_ALL, H5S_ALL, H5P_DEFAULT,
                          data.data()));
    }
    template <class T>
    void Read(std::vector<T> &data, DataSpaceSharedPtr filespace,
              PListSharedPtr dxpl = PList::Default())
    {
        DataTypeSharedPtr mem_t = DataTypeTraits<T>::GetType();
        hsize_t len;
        len = H5Sget_select_npoints(filespace->GetId());

        data.resize(len);

        DataSpaceSharedPtr memspace = DataSpace::OneD(len);
        H5_CALL(H5Dread, (m_Id, mem_t->GetId(), memspace->GetId(),
                          filespace->GetId(), dxpl->GetId(), data.data()));
    }
    template <class T>
    void Read(std::vector<T> &data, DataSpaceSharedPtr filespace,
              std::vector<std::vector<int>> &coords,
              PListSharedPtr dxpl = PList::Default())
    {
        DataTypeSharedPtr mem_t = DataTypeTraits<T>::GetType();
        hsize_t len;

        int w = coords[0].size();

        hsize_t *cds = new hsize_t[coords.size() * w];
        for (int i = 0; i < coords.size(); i++)
        {
            for (int j = 0; j < coords[i].size(); j++)
            {
                cds[i * w + j] = hsize_t(coords[i][j]);
            }
        }

        H5Sselect_elements(filespace->GetId(), H5S_SELECT_SET, coords.size(),
                           cds);

        delete[] cds;

        len                         = H5Sget_select_npoints(filespace->GetId());
        DataSpaceSharedPtr memspace = DataSpace::OneD(len);

        data.resize(len);

        H5_CALL(H5Dread, (m_Id, mem_t->GetId(), memspace->GetId(),
                          filespace->GetId(), dxpl->GetId(), data.data()));

        H5Sselect_all(filespace->GetId());
    }

    void ReadVectorString(std::vector<std::string> &data,
                          DataSpaceSharedPtr filespace,
                          PListSharedPtr dxpl = PList::Default())
    {
        char **rdata;
        DataTypeSharedPtr tp      = DataType::String();
        std::vector<hsize_t> dims = filespace->GetDims();
        rdata                     = (char **)malloc(dims[0] * sizeof(char *));

        H5_CALL(H5Dread, (m_Id, tp->GetId(), H5S_ALL, filespace->GetId(),
                          dxpl->GetId(), rdata));

        for (int i = 0; i < dims[0]; i++)
        {
            data.push_back(std::string(rdata[i]));
        }
        free(rdata);
    }

protected:
    LIB_UTILITIES_EXPORT void v_Close() override;

private:
    DataSet(hid_t id);
    friend class CanHaveGroupsDataSets;
};

template <class T>
typename DataTypeConversionPolicy<T>::ConvertedType DataTypeConversionPolicy<
    T>::Convert(const T &obj)
{
    return obj;
}

template <class T>
T DataTypeConversionPolicy<T>::Deconvert(
    const typename DataTypeConversionPolicy<T>::ConvertedType &obj)
{
    return obj;
}

template <class T> DataTypeSharedPtr DataTypeTraits<T>::GetType()
{
    return PredefinedDataType::Native<T>();
}

template <class T>
const void *DataTypeTraits<T>::GetAddress(
    const DataTypeTraits<T>::ConvertedType &obj)
{
    return &obj;
}

template <class T>
void *DataTypeTraits<T>::GetAddress(DataTypeTraits<T>::ConvertedType &obj)
{
    return &obj;
}

template <class T>
typename DataTypeTraits<T>::ConvertedType DataTypeTraits<T>::Convert(
    const T &obj)
{
    return Converter::Convert(obj);
}

template <class T>
T DataTypeTraits<T>::Deconvert(
    const typename DataTypeTraits<T>::ConvertedType &obj)
{
    return Converter::Deconvert(obj);
}
template <> struct DataTypeConversionPolicy<std::string>
{
    static const bool MustConvert = true;
    typedef const char *ConvertedType;
    typedef const char *ConvertedVectorElemType;
    inline static ConvertedType Convert(const std::string &obj)
    {
        return obj.c_str();
    }
    inline static std::string Deconvert(const ConvertedType &obj)
    {
        return std::string(obj);
    }
};

template <> inline DataTypeSharedPtr DataTypeTraits<std::string>::GetType()
{
    return DataType::String();
}

template <class T> DataTypeSharedPtr PredefinedDataType::Native()
{
    return DataTypeSharedPtr(
        new PredefinedDataType(DataTypeTraits<T>::NativeType));
}

template <typename T>
void CanHaveAttributes::SetAttribute(const std::string &name, const T &value)
{
    DataTypeSharedPtr type   = DataTypeTraits<T>::GetType();
    DataSpaceSharedPtr space = DataSpace::Scalar();
    AttributeSharedPtr attr  = CreateAttribute(name, type, space);

    typename DataTypeTraits<T>::ConvertedType conv =
        DataTypeTraits<T>::Convert(value);
    H5_CALL(H5Awrite, (attr->GetId(), type->GetId(),
                       DataTypeTraits<T>::GetAddress(conv)));
}

template <typename T>
void CanHaveAttributes::SetAttribute(const std::string &name,
                                     const std::vector<T> &value)
{
    typedef std::vector<
        typename DataTypeConversionPolicy<T>::ConvertedVectorElemType>
        Vec;
    Vec converted_vals;
    DataTypeSharedPtr type   = DataTypeTraits<T>::GetType();
    DataSpaceSharedPtr space = DataSpace::OneD(value.size());
    AttributeSharedPtr attr  = CreateAttribute(name, type, space);

    const void *converted_buf = nullptr;
    if (DataTypeConversionPolicy<T>::MustConvert)
    {
        converted_vals.resize(value.size());
        for (size_t i = 0; i < value.size(); ++i)
        {
            converted_vals[i] = DataTypeConversionPolicy<T>::Convert(value[i]);
        }
        converted_buf = &converted_vals[0];
    }
    else
    {
        converted_buf = &value[0];
    }

    H5_CALL(H5Awrite, (attr->GetId(), type->GetId(), converted_buf));
}

template <typename T>
void CanHaveAttributes::GetAttribute(const std::string &name, T &value)
{
    DataTypeSharedPtr type   = DataTypeTraits<T>::GetType();
    DataSpaceSharedPtr space = DataSpace::Scalar();
    AttributeSharedPtr attr  = OpenAttribute(name);

    typename DataTypeTraits<T>::ConvertedType conv;

    H5_CALL(H5Aread, (attr->GetId(), type->GetId(),
                      DataTypeTraits<T>::GetAddress(conv)));
    value = DataTypeTraits<T>::Deconvert(conv);
}

template <typename T>
void CanHaveAttributes::GetAttribute(const std::string &name,
                                     std::vector<T> &value)
{
    typedef std::vector<
        typename DataTypeConversionPolicy<T>::ConvertedVectorElemType>
        Vec;
    Vec converted_vals;
    DataTypeSharedPtr type   = DataTypeTraits<T>::GetType();
    AttributeSharedPtr attr  = OpenAttribute(name);
    DataSpaceSharedPtr space = attr->GetSpace();
    ASSERTL0(H5Sget_simple_extent_ndims(space->GetId()) == 1,
             "vector data not 1D");
    hsize_t len, maxdim;
    H5Sget_simple_extent_dims(space->GetId(), &len, &maxdim);

    value.resize(len);
    void *converted_buf = nullptr;
    if (DataTypeConversionPolicy<T>::MustConvert)
    {
        converted_vals.resize(len);
        converted_buf = &converted_vals[0];
    }
    else
    {
        converted_buf = &value[0];
    }

    H5_CALL(H5Aread, (attr->GetId(), type->GetId(), converted_buf));

    if (DataTypeConversionPolicy<T>::MustConvert)
    {
        typename Vec::iterator src             = converted_vals.begin(),
                               end             = converted_vals.end();
        typename std::vector<T>::iterator dest = value.begin();
        for (; src != end; ++src, ++dest)
        {
            *dest = DataTypeTraits<T>::Deconvert(*src);
        }
    }
}

template <class T>
DataSetSharedPtr CanHaveGroupsDataSets::CreateWriteDataSet(
    const std::string &name, const std::vector<T> &data,
    PListSharedPtr createPL, PListSharedPtr accessPL)
{
    DataTypeSharedPtr type   = DataTypeTraits<T>::GetType();
    DataSpaceSharedPtr space = DataSpace::OneD(data.size());
    DataSetSharedPtr dataset =
        CreateDataSet(name, type, space, createPL, accessPL);
    dataset->Write(data);
    return dataset;
}
} // namespace Nektar::LibUtilities::H5

#endif
