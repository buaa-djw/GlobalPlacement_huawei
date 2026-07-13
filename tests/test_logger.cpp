#include "utils/Logger.h"

#include <cassert>
#include <filesystem>
#include <fstream>
#include <iostream>
#include <sstream>
#include <string>

std::string readFile(const std::filesystem::path& p){ std::ifstream in(p); return std::string((std::istreambuf_iterator<char>(in)),{}); }

int main(){
    const auto dir = std::filesystem::temp_directory_path()/"gp_logger_test"; std::filesystem::remove_all(dir);
    const auto log = dir/"nested"/"test.log";
    Logger::instance().initialize(log.string(), LogLevel::Debug, true); assert(std::filesystem::exists(log));
    Logger::instance().log(LogLevel::Info,"hello info","unit.cpp",42); Logger::instance().flush(); auto s=readFile(log); assert(s.find("INFO")!=std::string::npos && s.find("hello info")!=std::string::npos && s.find("unit.cpp:42")!=std::string::npos);
    Logger::instance().initialize(log.string(), LogLevel::Info, false); LOG_DEBUG("hidden debug"); LOG_INFO("visible info"); Logger::instance().flush(); s=readFile(log); assert(s.find("hidden debug")==std::string::npos && s.find("visible info")!=std::string::npos);
    LOG_ERROR("error text"); Logger::instance().flush(); s=readFile(log); assert(s.find("ERROR")!=std::string::npos && s.find("error text")!=std::string::npos);
    bool fatal=false; try{ LOG_FATAL("fatal text"); }catch(const FatalLogError&){ fatal=true; } assert(fatal); s=readFile(log); assert(s.find("FATAL")!=std::string::npos && s.find("fatal text")!=std::string::npos);
    std::ostringstream cout_capture, cerr_capture; auto* old_out=std::cout.rdbuf(cout_capture.rdbuf()); auto* old_err=std::cerr.rdbuf(cerr_capture.rdbuf());
    Logger::instance().initialize(log.string(), LogLevel::Debug, true); LOG_INFO("console info"); LOG_ERROR("console error"); Logger::instance().flush(); std::cout.rdbuf(old_out); std::cerr.rdbuf(old_err);
    s=readFile(log); assert(s.find("console info")!=std::string::npos && s.find("console error")!=std::string::npos); assert(cout_capture.str().find("console info")!=std::string::npos); assert(cerr_capture.str().find("console error")!=std::string::npos);
    Logger::instance().initialize(log.string(), LogLevel::Info, false); Logger::instance().initialize(log.string(), LogLevel::Info, false); Logger::instance().shutdown(); Logger::instance().flush();
    return 0;
}
