#pragma once

#include <fstream>
#include <mutex>
#include <sstream>
#include <stdexcept>
#include <string>

/** @brief Severity threshold for Logger output. */
enum class LogLevel { Debug = 0, Info, Warn, Error, Fatal };

/** @brief Exception thrown after a FATAL log has already been emitted and flushed. */
class FatalLogError : public std::runtime_error { public: using std::runtime_error::runtime_error; };

/**
 * @brief Thread-safe singleton logger writing identical lines to console and file.
 *
 * A mutex protects initialization, shutdown, flushing, and both output sinks.
 * fatal() logs once, flushes while the state is consistent, then throws instead
 * of aborting so tests and main-level cleanup can run deterministically.
 */
class Logger {
public:
    static Logger& instance();
    void initialize(const std::string& file_path, LogLevel minimum_level = LogLevel::Info, bool enable_console = true);
    void shutdown();
    void flush();
    void log(LogLevel level, const std::string& message, const char* file, int line);
    [[noreturn]] void fatal(const std::string& message, const char* file, int line);
    bool initialized() const;
private:
    Logger() = default;
    std::string formatLine(LogLevel level, const std::string& message, const char* file, int line) const;
    mutable std::mutex mutex_;
    std::ofstream file_;
    LogLevel minimum_level_ = LogLevel::Info;
    bool enable_console_ = true;
    bool initialized_ = false;
};

LogLevel parseLogLevel(const std::string& text);
const char* logLevelName(LogLevel level);

#define LOG_STREAM_IMPL(level, expr) do { std::ostringstream _gp_log_os; _gp_log_os << expr; Logger::instance().log((level), _gp_log_os.str(), __FILE__, __LINE__); } while (false)
#define LOG_DEBUG(expr) LOG_STREAM_IMPL(LogLevel::Debug, expr)
#define LOG_INFO(expr) LOG_STREAM_IMPL(LogLevel::Info, expr)
#define LOG_WARN(expr) LOG_STREAM_IMPL(LogLevel::Warn, expr)
#define LOG_ERROR(expr) LOG_STREAM_IMPL(LogLevel::Error, expr)
#define LOG_FATAL(expr) do { std::ostringstream _gp_log_os; _gp_log_os << expr; Logger::instance().fatal(_gp_log_os.str(), __FILE__, __LINE__); } while (false)
