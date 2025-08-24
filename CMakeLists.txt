cmake_minimum_required(VERSION 3.29)
project(OOP_Phase2_complete)

set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

add_compile_definitions(GLEW_STATIC)

add_executable(OOP_Phase2_complete
        # Your main project file
        ProjectMainFile.cpp

        # Dear ImGui source files
        dependencies/imgui-master/imgui.cpp
        dependencies/imgui-master/imgui_draw.cpp
        dependencies/imgui-master/imgui_tables.cpp
        dependencies/imgui-master/imgui_widgets.cpp

        # NEW: Add the demo window source file
        dependencies/imgui-master/imgui_demo.cpp

        # Backend source files
        dependencies/imgui-master/backends/imgui_impl_glfw.cpp
        dependencies/imgui-master/backends/imgui_impl_opengl3.cpp

        # ImPlot source files
        dependencies/implot-master/implot.cpp
        dependencies/implot-master/implot_items.cpp
        save_load.h
)

target_include_directories(OOP_Phase2_complete PRIVATE
        # NEW: Add the vendor directory for Cereal and PFD
        ${CMAKE_CURRENT_SOURCE_DIR}/vendor

        # Existing directories
        ${CMAKE_CURRENT_SOURCE_DIR}/dependencies/imgui-master
        ${CMAKE_CURRENT_SOURCE_DIR}/dependencies/imgui-master/backends
        ${CMAKE_CURRENT_SOURCE_DIR}/dependencies/implot-master
        ${CMAKE_CURRENT_SOURCE_DIR}/dependencies/glfw-3.4.bin.WIN64/include
        ${CMAKE_CURRENT_SOURCE_DIR}/dependencies/glew-2.1.0/include
)

target_link_directories(OOP_Phase2_complete PRIVATE
        ${CMAKE_CURRENT_SOURCE_DIR}/dependencies/glfw-3.4.bin.WIN64/lib-mingw-w64
        ${CMAKE_CURRENT_SOURCE_DIR}/dependencies/glew-2.1.0/lib/Release/x64
)

target_link_libraries(OOP_Phase2_complete PRIVATE
        glew32s
        glfw3
        gdi32
        opengl32
)

# NEW: Link system libraries needed for Portable File Dialogs on Windows
if(WIN32)
    target_link_libraries(OOP_Phase2_complete PRIVATE
            Comdlg32
            Ole32
            Shell32
    )
endif()