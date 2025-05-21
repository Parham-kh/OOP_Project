cmake_minimum_required(VERSION 3.29)
project(OOP_Project)

set(CMAKE_CXX_STANDARD 14)

include_directories(.)

add_executable(OOP_Project
    main.cpp
        Exceptions/Exceptions.cpp
        Model/Node.cpp
        Model/Element.cpp
        Model/Resistor.cpp
        ProgramController/Circuit.cpp
        View/View.cpp
        Model/VoltageSource.cpp
        Model/Inductor.cpp
        Model/Capacitor.cpp
        Model/CurrentSource.cpp
)
