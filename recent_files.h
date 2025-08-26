#pragma once
#include <string>
#include <vector>
#include <algorithm>
#include <cstdlib>   // getenv
#include <cctype>
#include <fstream>
#include <filesystem>
#include "save_load.h"   // cereal already included in your project
// همچنین make_dirs(...) که قبلاً داری

// حداکثر تعداد آیتم‌ها
static const size_t kMaxRecent = 10;

inline void ensure_directory_exists(const std::string& path_str) {
    std::filesystem::path p(path_str);
    if (p.has_parent_path()) {
        std::filesystem::create_directories(p.parent_path());
    }
}

// مسیر ذخیره‌سازی: %APPDATA%\Opulator\recent.json  (اگر APPDATA نبود، پوشه‌ی کاری)
inline std::string recent_store_path() {
#ifdef _WIN32
    const char* appdata = std::getenv("APPDATA");
    if (appdata && *appdata) {
        std::string dir = std::string(appdata) + "\\Opulator";
        return dir + "\\recent.json";
    }
#endif
    return "recent/recent.json";
}

// نرمال‌سازی سبک مسیر برای مقایسه (حذف اسپیس‌های دو سر، یکسان‌سازی \ / ، و lowercase)
inline std::string normalize_path_for_compare(std::string p) {
    auto trim = [](std::string& s){
        while (!s.empty() && (s.back()==' '||s.back()=='\t')) s.pop_back();
        size_t i=0; while (i<s.size() && (s[i]==' '||s[i]=='\t')) ++i;
        if (i) s.erase(0,i);
    };
    trim(p);
    for (char& ch : p) {
        if (ch == '/') ch = '\\';
        ch = (char)std::tolower((unsigned char)ch);
    }
    return p;
}

// ساختار ذخیره‌سازی برای cereal
struct RecentStore {
    int version{1};
    std::vector<std::string> files; // UTF-8 full paths

    template<class Ar>
    void serialize(Ar& ar) {
        ar( CEREAL_NVP(version), CEREAL_NVP(files) );
    }
};

inline void recent_load(std::vector<std::string>& out) {
    out.clear();
    std::ifstream is(recent_store_path());
    if (!is) return;
    try {
        cereal::JSONInputArchive ar(is);
        RecentStore rs; ar(cereal::make_nvp("recent", rs));
        out = rs.files;
    } catch (...) { /* ignore corrupt */ }
}

inline void recent_save(const std::vector<std::string>& in) {
    std::string path = recent_store_path();
    ensure_directory_exists(path); // Create the directory if it doesn't exist

    RecentStore rs; rs.files = in;
    std::ofstream os(path);
    if (!os) return;
    cereal::JSONOutputArchive ar(os);
    ar(cereal::make_nvp("recent", rs));
}

// افزودن آیتم جدید به لیست (حذف دابلکِیت و محدود کردن سایز)
inline void recent_add(const std::string& path) {
    std::vector<std::string> cur; recent_load(cur);
    const std::string key = normalize_path_for_compare(path);

    // حذف دابلکِیت (case-insensitive + slash fix)
    cur.erase(std::remove_if(cur.begin(), cur.end(), [&](const std::string& s){
        return normalize_path_for_compare(s) == key;
    }), cur.end());

    // افزودن به ابتدای لیست
    cur.insert(cur.begin(), path);

    // محدود کردن سایز
    if (cur.size() > kMaxRecent) cur.resize(kMaxRecent);

    recent_save(cur);
}

inline std::vector<std::string> recent_list() {
    std::vector<std::string> cur; recent_load(cur);
    return cur;
}

inline void recent_clear() {
    std::vector<std::string> empty;
    recent_save(empty);
}
