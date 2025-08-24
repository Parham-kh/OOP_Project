#pragma once
#include <vector>
#include <string>

// cereal
#include <cereal/archives/json.hpp>
#include <cereal/types/vector.hpp>
#include <cereal/types/string.hpp>

struct SNode {
    int id{}; std::string name; int x{0}, y{0};
    template<class Ar> void serialize(Ar& ar){
        ar(CEREAL_NVP(id), CEREAL_NVP(name), CEREAL_NVP(x), CEREAL_NVP(y));
    }
};

struct Point {
    int x{0}, y{0};
    template<class Ar> void serialize(Ar& ar) {
        ar(CEREAL_NVP(x), CEREAL_NVP(y));
    }
};

struct SWire {
    int id{}; std::vector<Point> points; // polyline
    template<class Ar> void serialize(Ar& ar){ ar(CEREAL_NVP(id), CEREAL_NVP(points)); }
};

struct SElement {
    std::string kind, name;
    int n1{-1}, n2{-1};
    double value{0.0};
    int x{0}, y{0};
    bool isVertical{false};

    // Source parameters
    std::string sourceType{"DC"}; // "DC", "SINE", "PULSE", "DELTA", "VCVS", "VCCS" etc.
    double amp{0.0}, freq{0.0};

    // Pulse parameters
    double v_initial{0.0}, v_on{0.0};
    double i_initial{0.0}, i_on{0.0};
    double t_delay{0.0}, t_rise{0.0}, t_fall{0.0}, t_on{0.0}, t_period{0.0};

    // Delta parameters
    double t_pulse{0.0}, area{1.0};

    // Dependent source parameters
    double gain{1.0};
    std::string ctrlNodeP_name, ctrlNodeN_name;

    template<class Ar> void serialize(Ar& ar){
        ar(CEREAL_NVP(kind), CEREAL_NVP(name), CEREAL_NVP(n1), CEREAL_NVP(n2),
           CEREAL_NVP(value), CEREAL_NVP(x), CEREAL_NVP(y), CEREAL_NVP(isVertical));

        // Conditionally save source parameters
        if (kind == "V" || kind == "I") {
            ar(CEREAL_NVP(sourceType));
            if (sourceType == "SINE") { ar(CEREAL_NVP(amp), CEREAL_NVP(freq)); }
            else if (sourceType == "PULSE") {
                ar(CEREAL_NVP(v_initial), CEREAL_NVP(v_on), CEREAL_NVP(i_initial), CEREAL_NVP(i_on),
                   CEREAL_NVP(t_delay), CEREAL_NVP(t_rise), CEREAL_NVP(t_fall), CEREAL_NVP(t_on), CEREAL_NVP(t_period));
            }
            else if (sourceType == "DELTA") { ar(CEREAL_NVP(t_pulse), CEREAL_NVP(area)); }
            else if (sourceType == "VCVS" || sourceType == "VCCS") {
                ar(CEREAL_NVP(gain), CEREAL_NVP(ctrlNodeP_name), CEREAL_NVP(ctrlNodeN_name));
            }
        }
    }
};

struct TransientCfg {
    double t_stop{1e-3}, t0{0}, t_step{1e-6};
    template<class Ar> void serialize(Ar& ar){ ar(CEREAL_NVP(t_stop), CEREAL_NVP(t0), CEREAL_NVP(t_step)); }
};

struct ACSweepCfg {
    std::string type{"Linear"}; double w_start{10}, w_stop{1e6}; int points{200};
    template<class Ar> void serialize(Ar& ar){
        ar(CEREAL_NVP(type), CEREAL_NVP(w_start), CEREAL_NVP(w_stop), CEREAL_NVP(points));
    }
};

struct PhaseSweepCfg {
    double w_base{1000}, phi_start{0}, phi_stop{3.14159}; int points{181};
    template<class Ar> void serialize(Ar& ar){
        ar(CEREAL_NVP(w_base), CEREAL_NVP(phi_start), CEREAL_NVP(phi_stop), CEREAL_NVP(points));
    }
};

struct ProjectSettings {
    TransientCfg transient; ACSweepCfg ac; PhaseSweepCfg phase;
    template<class Ar> void serialize(Ar& ar){
        ar(CEREAL_NVP(transient), CEREAL_NVP(ac), CEREAL_NVP(phase));
    }
};

struct SCircuit {
    int file_version{1};
    std::vector<SNode>    Snodes;
    std::vector<SWire>    Swires;
    std::vector<SElement> Selements;
    ProjectSettings       settings;

    template<class Ar> void serialize(Ar& ar){
        ar( cereal::make_nvp("version",  file_version),
            cereal::make_nvp("nodes",    Snodes),
            cereal::make_nvp("wires",    Swires),
            cereal::make_nvp("elements", Selements),
            cereal::make_nvp("settings", settings) );
    }
};

bool saveProject(const std::string& folderPath, const SCircuit& ckt);
bool loadProject(const std::string& folderPath, SCircuit& out);
