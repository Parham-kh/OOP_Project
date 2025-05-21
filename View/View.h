#ifndef VIEW_H
#define VIEW_H
#pragma once
#include <string>
#include <regex>
#include <iostream>
#include "ProgramController/Circuit.h"
using namespace std;
using namespace Controller;

namespace View {
    class view {
        Circuit* circuit;

    public:
        void run();
    };
}

#endif //VIEW_H
