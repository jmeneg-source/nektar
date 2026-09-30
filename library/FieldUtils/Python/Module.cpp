///////////////////////////////////////////////////////////////////////////////
//
// File: Module.cpp
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
// Description: Python wrapper for Module.
//
///////////////////////////////////////////////////////////////////////////////

#include <FieldUtils/Module.h>
#include <LibUtilities/Python/BasicUtils/NekFactory.hpp>
#include <LibUtilities/Python/BasicUtils/SharedArray.hpp>
#include <LibUtilities/Python/NekPyConfig.hpp>

#include <boost/program_options.hpp>

#ifdef NEKTAR_USING_VTK
#include <FieldUtils/OutputModules/OutputVtk.h>
#include <FieldUtils/Python/VtkWrapper.hpp>
#endif

using namespace Nektar;
using namespace Nektar::FieldUtils;

/**
 * @brief Module wrapper to handle virtual function calls in @c Module and its
 * subclasses as defined by the template parameter @tparam MODTYPE.
 */
#pragma GCC visibility push(hidden)
template <typename MODTYPE>
struct ModuleWrap : public MODTYPE, public py::trampoline_self_life_support
{
    /**
     * @brief Constructor, which is identical to FieldUtils::Module.
     *
     * @param field  Input field.
     */
    ModuleWrap(FieldSharedPtr field) : MODTYPE(field)
    {
    }

    /**
     * @brief Concrete implementation of the Module::Process function.
     */
    void v_Process([[maybe_unused]] po::variables_map &vm) override
    {
        py::gil_scoped_acquire gil;
        py::function override = py::get_override(this, "Process");

        if (this->m_f->m_nParts > 1)
        {
            if (this->GetModulePriority() == eOutput)
            {
                this->m_f->m_comm = this->m_f->m_partComm;
                if (this->GetModuleName() != "OutputInfo")
                {
                    this->RegisterConfig("writemultiplefiles");
                }
            }
            else if (this->GetModulePriority() == eCreateGraph)
            {
                this->m_f->m_comm = this->m_f->m_partComm;
            }
            else
            {
                this->m_f->m_comm = this->m_f->m_defComm;
            }
        }

        if (override)
        {
            override();
            return;
        }

        throw ErrorUtil::NekError("Process() is a pure virtual function");
    }

    std::string v_GetModuleName() override
    {
        py::gil_scoped_acquire gil;
        py::function override = py::get_override(this, "GetModuleName");

        if (override)
        {
            return py::cast<std::string>(override());
        }

        throw ErrorUtil::NekError("GetModuleName() is a pure virtual function");
    }

    ModulePriority v_GetModulePriority() override
    {
        py::gil_scoped_acquire gil;
        py::function override = py::get_override(this, "GetModulePriority");

        if (override)
        {
            return py::cast<ModulePriority>(override());
        }

        throw ErrorUtil::NekError(
            "GetModulePriority() is a pure virtual function");
    }

    /**
     * @brief Defines a configuration option for this module.
     *
     * @param key     The name of the configuration option.
     * @param def     The option's default value.
     * @param desc    A text description of the option.
     * @param isBool  If true, this option is a boolean-type (true/false).
     */
    void AddConfigOption(std::string key, std::string def, std::string desc,
                         bool isBool)
    {
        ConfigOption conf(isBool, def, desc);
        this->m_config[key] = conf;
    }

    // We expose Module::m_f as a public member variable so that we can
    // adjust this using Python attributes.
    using Module::m_f;
};
#pragma GCC visibility pop

template <typename MODTYPE> struct ModulePublic : public MODTYPE
{
    using MODTYPE::v_GetModuleName;
    using MODTYPE::v_GetModulePriority;
    using MODTYPE::v_Process;
};

// Wrapper around Module::Process(&vm).
// Performs switching of m_comm if nparts > 1.
void Module_Process(ModuleSharedPtr m)
{
    if (m->m_f->m_nParts > 1)
    {
        if (m->GetModulePriority() == eOutput)
        {
            m->m_f->m_comm = m->m_f->m_partComm;
            if (m->GetModuleName() != "OutputInfo")
            {
                m->RegisterConfig("writemultiplefiles");
            }
        }
        else if (m->GetModulePriority() == eCreateGraph)
        {
            m->m_f->m_comm = m->m_f->m_partComm;
        }
        else
        {
            m->m_f->m_comm = m->m_f->m_defComm;
        }
    }
    m->SetDefaults();
    m->Process(m->m_f->m_vm);
}

template <typename T>
T Module_GetConfig(std::shared_ptr<Module> mod, const std::string &key)
{
    return mod->GetConfigOption(key).as<T>();
}

template <typename MODTYPE> struct ModuleTypeProxy
{
};

template <> struct ModuleTypeProxy<InputModule>
{
    static const ModuleType value = eInputModule;
};

template <> struct ModuleTypeProxy<ProcessModule>
{
    static const ModuleType value = eProcessModule;
};

template <> struct ModuleTypeProxy<OutputModule>
{
    static const ModuleType value = eOutputModule;
};

const ModuleType ModuleTypeProxy<InputModule>::value;
const ModuleType ModuleTypeProxy<ProcessModule>::value;
const ModuleType ModuleTypeProxy<OutputModule>::value;

#ifdef NEKTAR_USING_VTK
vtkUnstructuredGrid *Module_GetVtkGrid(std::shared_ptr<Module> mod)
{
    std::shared_ptr<OutputVtk> vtkModule =
        std::dynamic_pointer_cast<OutputVtk>(mod);

    using NekError = ErrorUtil::NekError;

    if (!vtkModule)
    {
        throw NekError("This module is not the OutputVtk module, cannot get"
                       "VTK data.");
    }

    return vtkModule->GetVtkGrid();
}
#else
void Module_GetVtkGrid(std::shared_ptr<Module> mod)
{
    using NekError = ErrorUtil::NekError;
    throw NekError("Nektar++ has not been compiled with VTK support.");
}
#endif

/**
 * @brief Lightweight wrapper for Module factory creation function.
 *
 * @param  modType  Module type (input/process/output).
 * @param  modName  Module name (typically filename extension).
 * @param  field    Field that will be passed between modules.
 * @tparam MODTYPE  Subclass of Module (e.g #InputModule, #OutputModule)
 */
template <typename MODTYPE>
ModuleSharedPtr Module_Create(py::args args, const py::kwargs &kwargs)
{
    ModuleType modType = ModuleTypeProxy<MODTYPE>::value;

    using NekError = ErrorUtil::NekError;

    if (modType == eProcessModule && py::len(args) != 2)
    {
        throw NekError("ProcessModule.Create() requires two arguments: "
                       "module name and a Field object.");
    }
    else if (modType != eProcessModule && py::len(args) < 2)
    {
        throw NekError(ModuleTypeMap[modType] +
                       "Module.Create() requires "
                       "two arguments: module name and a Field object; "
                       "optionally a filename.");
    }

    std::string modName = py::str(args[0]);
    ModuleKey modKey    = std::make_pair(modType, modName);

    FieldSharedPtr field;

    try
    {
        field = py::cast<FieldSharedPtr>(args[1]);
    }
    catch (...)
    {
        throw NekError("Second argument to Create() should be a Field object.");
    }

    if (!GetModuleFactory().ModuleExists(modKey))
    {
        throw ErrorUtil::NekError("Module '" + modName + "' does not exist.");
    }

    ModuleSharedPtr mod = GetModuleFactory().CreateInstance(modKey, field);

    if (modType == eInputModule)
    {
        // For input modules we can try to interpret the remaining arguments as
        // input files. Assume that the file's type is identical to the module
        // name.
        for (int i = 2; i < py::len(args); ++i)
        {
            std::string in_fname = py::str(args[i]);
            mod->RegisterConfig("infile", in_fname);
            mod->AddFile(modName, in_fname);
        }
    }
    else if (modType == eOutputModule && py::len(args) >= 3)
    {
        // For output modules we can try to interpret the remaining argument as
        // an output file.
        mod->RegisterConfig("outfile", py::str(args[2]));
    }

    // Process keyword arguments.
    for (auto &kv : kwargs)
    {
        std::string arg = py::str(kv.first);

        if (arg == "infile" && modKey.first == eInputModule)
        {
            if (!py::isinstance<py::dict>(kv.second))
            {
                throw NekError("infile should be a dictionary.");
            }

            py::dict ftype_fname_dict = py::cast<py::dict>(kv.second);

            for (auto &kv2 : ftype_fname_dict)
            {
                std::string f_type = py::str(kv2.first);
                std::string f_name = py::str(kv2.second);
                mod->RegisterConfig(arg, f_name);
                mod->AddFile(f_type, f_name);
            }
        }
        else
        {
            std::string val = py::str(kv.second);
            mod->RegisterConfig(arg, val);
        }
    }

    mod->SetDefaults();

    return mod;
}

/**
 * @brief Lightweight wrapper for FieldUtils::Module::RegisterConfig.
 *
 * @param mod    Module to call
 * @param key    Configuration key.
 * @param value  Optional value (some configuration options are boolean).
 */
void Module_RegisterConfig(std::shared_ptr<Module> mod, std::string const &key,
                           std::string const &value)
{
    mod->RegisterConfig(key, value);
}

template <typename MODTYPE>
void ModuleWrap_AddConfigOption(std::shared_ptr<ModuleWrap<MODTYPE>> mod,
                                std::string const &key,
                                std::string const &defValue,
                                std::string const &desc, bool isBool)
{
    mod->AddConfigOption(key, defValue, desc, isBool);
}

template <typename MODTYPE> struct PythonModuleClass
{
    PythonModuleClass(py::module &m, std::string modName)
    {
        py::classh<MODTYPE, Module, ModuleWrap<MODTYPE>>(m, modName.c_str())
            .def(py::init<FieldSharedPtr>())

            .def("AddConfigOption", &ModuleWrap_AddConfigOption<MODTYPE>,
                 py::arg("key"), py::arg("defValue"), py::arg("desc"),
                 py::arg("isBool") = false)

            // Allow direct access to field object through a property.
            .def_readwrite("field", &ModuleWrap<MODTYPE>::m_f)

            // Process function for this module.
            .def("Process", &ModulePublic<MODTYPE>::v_Process)
            .def("Run", &ModulePublic<MODTYPE>::v_Process)
            .def_static("Create", &Module_Create<MODTYPE>);
    }
};

void export_Module(py::module &m)
{
    static NekFactory_Register<ModuleFactory> fac(GetModuleFactory());

    // Export ModuleType enum.
    NEKPY_WRAP_ENUM_STRING(m, ModuleType, ModuleTypeMap);

    // Wrapper for the Module class. Note that since Module contains a pure
    // virtual function, we need the ModuleWrap helper class to handle this for
    // us. In the lightweight wrappers above, we therefore need to ensure we're
    // passing std::shared_ptr<Module> as the first argument, otherwise they
    // won't accept objects constructed from Python.
    py::classh<Module, ModuleWrap<Module>>(m, "Module")
        .def(py::init<FieldSharedPtr>())

        // Process function for this module.
        .def("Process", &Module_Process)
        .def("Run", &Module_Process)

        // Configuration options.
        .def("RegisterConfig", Module_RegisterConfig, py::arg("key"),
             py::arg("value") = "")
        .def("PrintConfig", &Module::PrintConfig)
        .def("SetDefaults", &Module::SetDefaults)
        .def("GetStringConfig", Module_GetConfig<std::string>)
        .def("GetFloatConfig", Module_GetConfig<double>)
        .def("GetIntConfig", Module_GetConfig<int>)
        .def("GetBoolConfig", Module_GetConfig<bool>)
        .def("AddConfigOption", &ModuleWrap_AddConfigOption<Module>,
             py::arg("key"), py::arg("defValue"), py::arg("desc"),
             py::arg("isBool") = false)
        .def("GetVtkGrid", &Module_GetVtkGrid)

        // Allow direct access to field object through a property.
        .def_readwrite("field", &ModuleWrap<Module>::m_f)

        // Factory functions.
        .def_static("Register", [](ModuleType const &modType,
                                   std::string const &modName,
                                   py::object &obj) {
            ModuleKey key(modType, modName);
            fac(key, obj, std::string(ModuleTypeMap[modType]) + "_" + modName);
        });

    PythonModuleClass<InputModule>(m, "InputModule");
    PythonModuleClass<ProcessModule>(m, "ProcessModule");
    PythonModuleClass<OutputModule>(m, "OutputModule");
}
