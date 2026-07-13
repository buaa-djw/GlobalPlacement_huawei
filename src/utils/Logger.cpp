#include "utils/Logger.h"

#include <chrono>
#include <ctime>
#include <filesystem>
#include <iomanip>
#include <iostream>

Logger& Logger::instance() { static Logger logger; return logger; }

const char* logLevelName(LogLevel level) {
    switch (level) { case LogLevel::Debug: return "DEBUG"; case LogLevel::Info: return "INFO"; case LogLevel::Warn: return "WARN"; case LogLevel::Error: return "ERROR"; case LogLevel::Fatal: return "FATAL"; }
    return "UNKNOWN";
}

LogLevel parseLogLevel(const std::string& text) {
    if (text == "debug") return LogLevel::Debug;
    if (text == "info") return LogLevel::Info;
    if (text == "warn") return LogLevel::Warn;
    if (text == "error") return LogLevel::Error;
    throw std::invalid_argument("unknown log level: " + text);
}

void Logger::initialize(const std::string& file_path, LogLevel minimum_level, bool enable_console) {
    std::lock_guard<std::mutex> lock(mutex_);
    if (file_.is_open()) file_.close();
    minimum_level_ = minimum_level;
    enable_console_ = enable_console;
    const std::filesystem::path path(file_path);
    if (path.has_parent_path()) std::filesystem::create_directories(path.parent_path());
    // Overwrite mode keeps each program run self-contained instead of mixing old diagnostics.
    file_.open(path, std::ios::out | std::ios::trunc);
    if (!file_) {
        std::cerr << "Cannot open log file: " << file_path << '\n';
        throw std::runtime_error("Logger::initialize: cannot open log file: " + file_path);
    }
    initialized_ = true;
}

void Logger::shutdown() {
    std::lock_guard<std::mutex> lock(mutex_);
    if (file_.is_open()) { file_.flush(); file_.close(); }
    initialized_ = false;
}

void Logger::flush() {
    std::lock_guard<std::mutex> lock(mutex_);
    if (file_.is_open()) file_.flush();
    std::cout.flush(); std::cerr.flush();
}

bool Logger::initialized() const { std::lock_guard<std::mutex> lock(mutex_); return initialized_; }

std::string Logger::formatLine(LogLevel level, const std::string& message, const char* file, int line) const {
    using clock = std::chrono::system_clock;
    const auto now = clock::now();
    const auto ms = std::chrono::duration_cast<std::chrono::milliseconds>(now.time_since_epoch()) % 1000;
    const std::time_t t = clock::to_time_t(now);
    std::tm tm{};
#if defined(_WIN32)
    localtime_s(&tm, &t);
#else
    localtime_r(&t, &tm);
#endif
    std::ostringstream os;
    os << '[' << std::put_time(&tm, "%Y-%m-%d %H:%M:%S") << '.' << std::setw(3) << std::setfill('0') << ms.count() << "] "
       << '[' << logLevelName(level) << "] [" << file << ':' << line << "] " << message;
    return os.str();
}

void Logger::log(LogLevel level, const std::string& message, const char* file, int line) {
    std::lock_guard<std::mutex> lock(mutex_);
    if (static_cast<int>(level) < static_cast<int>(minimum_level_)) return;
    const std::string out = formatLine(level, message, file, line);
    if (enable_console_) {
        (static_cast<int>(level) >= static_cast<int>(LogLevel::Error) ? std::cerr : std::cout) << out << '\n';
    }
    if (file_.is_open()) file_ << out << '\n';
}

[[noreturn]] void Logger::fatal(const std::string& message, const char* file, int line) {
    log(LogLevel::Fatal, message, file, line);
    flush();
    throw FatalLogError(message);
}
