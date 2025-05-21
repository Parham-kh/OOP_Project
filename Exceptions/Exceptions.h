
#ifndef EXCEPTIONS_H
#define EXCEPTIONS_H
#pragma once
#include <exception>

class ResistorInvalidAmount : public std::exception {
public:
    const char* what() const noexcept override {
        return "Error: Resistance cannot be zero or negative";
    }
};

class SyntaxError : public std::exception {
public:
    const char* what() const noexcept override {
        return "Error: Syntax error";
    }
};



#endif //EXCEPTIONS_H
