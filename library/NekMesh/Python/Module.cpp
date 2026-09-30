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

#include <LibUtilities/Python/BasicUtils/NekFactory.hpp>
#include <LibUtilities/Python/NekPyConfig.hpp>
#include <NekMesh/Module/Module.h>
#include <NekMesh/Python/NekMesh.h>

using namespace Nektar;
using namespace Nektar::NekMesh;

/**
 * @brief A basic Python logger object, which writes log messages to sys.stdout.
 *
 * This log streamer writes log messages to Python's sys.stdout. Although this
 * may seem superfluous, in some cases certain Python applications (in
 * particular, Jupyter notebooks), override sys.stdout to redirect output to
 * e.g. a browser or other IOStream. This therefore enables C++ log messages to
 * be redirected accordingly.
 */
class PythonStream : public LogOutput
{
public:
    /// Default constructor.
    PythonStream() : LogOutput()
    {
    }

    /**
     * @brief Write a log message @p msg.
     *
     * The logger runs the Python command `sys.stdout.write(msg)`.
     *
     * @param msg   The message to be printed.
     */
    void Log(const std::string &msg) override
    {
        // Create a Python string with the message.
        py::str pyMsg(msg);
        // Write to the Python stdout.
        py::module_::import("sys").attr("stdout").attr("write")(msg);
    }

private:
    /// Finalise function. Nothing required in this case.
    void Finalise() override
    {
    }
};

/// A horrible but necessary static variable to enable/disable verbose output by
/// default.
static bool default_verbose = false;

/**
 * @brief Module wrapper to handle virtual function calls in @c Module and its
 * subclasses as defined by the template parameter @tparam MODTYPE.
 */
#pragma GCC visibility push(hidden)
struct ModuleWrap : public Module, public py::trampoline_self_life_support
{
    /**
     * @brief Constructor, which is identical to NekMesh::Module::Module.
     *
     * @param mesh  Input mesh.
     */
    ModuleWrap(MeshSharedPtr mesh) : Module(mesh)
    {
    }

    /**
     * @brief Concrete implementation of the Module::Process function.
     */
    void Process() override
    {
        PYBIND11_OVERRIDE_PURE(void, Module, Process, );
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

    // We expose Module::m_mesh as a public member variable so that we can
    // adjust this using Python attributes.
    using Module::m_mesh;
};
#pragma GCC visibility pop

/**
 * @brief Module wrapper to handle virtual function calls in @c Module and its
 * subclasses as defined by the template parameter @tparam MODTYPE.
 */
#pragma GCC visibility push(hidden)
template <typename MODTYPE>
struct SubmoduleWrap : public MODTYPE, public py::trampoline_self_life_support
{
    /**
     * @brief Constructor, which is identical to NekMesh::Module::Module.
     *
     * @param mesh  Input mesh.
     */
    SubmoduleWrap(MeshSharedPtr mesh) : MODTYPE(mesh)
    {
    }

    /**
     * @brief Concrete implementation of the Module::Process function.
     */
    void Process() override
    {
        PYBIND11_OVERRIDE_PURE(void, MODTYPE, Process, );
    }
};
#pragma GCC visibility pop

template <typename T>
T Module_GetConfig(std::shared_ptr<Module> mod, const std::string &key)
{
    return mod->GetConfigOption(key).as<T>();
}

/**
 * @brief Lightweight wrapper for NekMesh::Module::ProcessEdges.
 *
 * @param mod        Module to call.
 * @param reprocess  If true then edges will be reprocessed (i.e. match 1D
 *                   elements to 2D element edges to construct boundaries).
 */
void Module_ProcessEdges(std::shared_ptr<Module> mod, bool reprocess = true)
{
    mod->ProcessEdges(reprocess);
}

/**
 * @brief Lightweight wrapper for NekMesh::Module::ProcessFaces.
 *
 * @param mod        Module to call.
 * @param reprocess  If true then faces will be reprocessed (i.e. match 2D
 *                   elements to 3D element faces to construct boundaries).
 */
void Module_ProcessFaces(std::shared_ptr<Module> mod, bool reprocess = true)
{
    mod->ProcessFaces(reprocess);
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

/**
 * @brief Lightweight wrapper for Module factory creation function.
 *
 * @param  modType  Module type (input/process/output).
 * @param  modName  Module name (typically filename extension).
 * @param  mesh     Mesh that will be passed between modules.
 * @tparam MODTYPE  Subclass of Module (e.g #InputModule, #OutputModule)
 */
template <typename MODTYPE>
ModuleSharedPtr Module_Create(py::args args, const py::kwargs &kwargs)
{
    ModuleType modType = ModuleTypeProxy<MODTYPE>::value;

    if (modType == eProcessModule && py::len(args) != 2)
    {
        throw NekMeshError("ProcessModule.Create() requires two arguments: "
                           "module name and a Mesh object.");
    }
    else if (modType != eProcessModule && py::len(args) != 3)
    {
        throw NekMeshError(ModuleTypeMap[modType] +
                           "Module.Create() requires "
                           "three arguments: module name, a Mesh object, and a "
                           "filename");
    }

    std::string modName = py::cast<std::string>(args[0]);
    ModuleKey modKey    = std::make_pair(modType, modName);

    if (!GetModuleFactory().ModuleExists(modKey))
    {
        throw NekMeshError("Module '" + modName + "' does not exist.");
    }

    MeshSharedPtr mesh;
    try
    {
        mesh = py::cast<MeshSharedPtr>(args[1]);
    }
    catch (...)
    {
        throw NekMeshError("Second argument to Create() should be a mesh "
                           "object.");
    }

    ModuleSharedPtr mod = GetModuleFactory().CreateInstance(modKey, mesh);

    // First argument for input/output module should be the filename.
    if (modKey.first == eInputModule)
    {
        mod->RegisterConfig("infile", py::cast<std::string>(args[2]));
    }
    else if (modKey.first == eOutputModule)
    {
        mod->RegisterConfig("outfile", py::cast<std::string>(args[2]));
    }

    // Set default verbosity.
    bool verbose = default_verbose;

    for (auto &kv : kwargs)
    {
        std::string arg = py::str(kv.first), val;

        // Enable or disable verbose for this module accordingly.
        if (arg == "verbose")
        {
            verbose = py::cast<bool>(kv.second);
            break;
        }
    }

    // Set a logger for this module.
    auto pythonLog = std::make_shared<PythonStream>();
    Logger log     = Logger(pythonLog, verbose ? VERBOSE : INFO);
    mod->SetLogger(log);

    for (auto &kv : kwargs)
    {
        std::string arg = py::str(kv.first), val;

        // Enable or disable verbose for this module accordingly.
        if (arg == "verbose")
        {
            continue;
        }

        val = py::str(kv.second);
        mod->RegisterConfig(arg, val);
    }

    // Set other default arguments.
    mod->SetDefaults();

    return mod;
}

/**
 * @brief Lightweight wrapper for NekMesh::Module::RegisterConfig.
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

/**
 * @brief Add a configuration option for the module.
 *
 * @tparam MODTYPE   Module type (e.g. #ProcessModule).
 * @param  mod       Module object.
 * @param  key       Name of the configuration option.
 * @param  defValue  Default value.
 * @param  desc      Description of the option.
 * @param  isBool    If true, denotes that the option will be bool-type.
 */
void ModuleWrap_AddConfigOption(std::shared_ptr<ModuleWrap> mod,
                                std::string const &key,
                                std::string const &defValue,
                                std::string const &desc, bool isBool)
{
    mod->AddConfigOption(key, defValue, desc, isBool);
}

/**
 * @brief Enables or disables verbose output by default.
 */
void Module_Verbose(bool verbose)
{
    default_verbose = verbose;
}

/**
 * @brief Wrapper for subclasses of the Module class, e.g. #InputModule, which
 * can then be inhereted from inside Python.
 */
template <typename MODTYPE> struct PythonModuleClass
{
    PythonModuleClass(py::module &m, std::string modName)
    {
        py::classh<MODTYPE, Module, SubmoduleWrap<MODTYPE>>(m, modName.c_str())
            .def(py::init<MeshSharedPtr>())
            .def("Create", &Module_Create<MODTYPE>);
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
    py::classh<Module, ModuleWrap>(m, "Module")
        .def(py::init<MeshSharedPtr>())

        // Process function for this module.
        .def("Process", &Module::Process)

        // Configuration options.
        .def("RegisterConfig", &Module_RegisterConfig, py::arg("key"),
             py::arg("value") = "")
        .def("PrintConfig", &Module::PrintConfig)
        .def("SetDefaults", &Module::SetDefaults)
        .def("GetStringConfig", Module_GetConfig<std::string>)
        .def("GetFloatConfig", Module_GetConfig<double>)
        .def("GetIntConfig", Module_GetConfig<int>)
        .def("GetBoolConfig", Module_GetConfig<bool>)
        .def("AddConfigOption", &ModuleWrap_AddConfigOption, py::arg("key"),
             py::arg("defValue"), py::arg("desc"), py::arg("isBool") = false)

        // Mesh accessor method.
        .def("GetMesh", &Module::GetMesh)

        // Mesh processing functions.
        .def("ProcessVertices", &Module::ProcessVertices)
        .def("ProcessEdges", &Module_ProcessEdges,
             py::arg("reprocessEdges") = true)
        .def("ProcessFaces", &Module_ProcessFaces,
             py::arg("reprocessFaces") = true)
        .def("ProcessElements", &Module::ProcessElements)
        .def("ProcessComposites", &Module::ProcessComposites)
        .def("ClearElementLinks", &Module::ClearElementLinks)

        // Allow direct access to mesh object through a property.
        .def_readwrite("mesh", &ModuleWrap::m_mesh)

        // Factory functions.
        .def_static("Register",
                    [](ModuleType const &modType, std::string const &modName,
                       py::object &obj) {
                        ModuleKey key(modType, modName);
                        fac(key, obj,
                            std::string(ModuleTypeMap[modType]) + "_" +
                                modName);
                    })

        // Enable verbose output (or not).
        .def_static("Verbose", &Module_Verbose);

    PythonModuleClass<InputModule>(m, "InputModule");
    PythonModuleClass<ProcessModule>(m, "ProcessModule");
    PythonModuleClass<OutputModule>(m, "OutputModule");
}
