#ifndef VIEW_H
#define VIEW_H
#pragma once
#include <string>
#include <regex>
#include <iostream>
#include "ProgramController/Circuit.h"

namespace View {
    class view {
        Controller::Circuit* circuit;

    public:
        view(Controller::Circuit* c);
        void run();
    };
}

#endif //VIEW_H
