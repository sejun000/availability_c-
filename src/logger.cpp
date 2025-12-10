#include "logger.hpp"
#include <chrono>
#include <ctime>
#include <iomanip>
#include <sstream>
#include <thread>

Logger::Logger() : min_level_(Level::INFO), console_enabled_(true) {}

Logger::~Logger() {
    if (log_file_.is_open()) {
        log_file_.close();
    }
}

Logger& Logger::getInstance() {
    static Logger instance;
    return instance;
}

void Logger::setLogFile(const std::string& filename) {
    std::lock_guard<std::mutex> lock(mutex_);
    if (log_file_.is_open()) {
        log_file_.close();
    }
    log_file_.open(filename, std::ios::app);
}

void Logger::setLevel(Level level) {
    std::lock_guard<std::mutex> lock(mutex_);
    min_level_ = level;
}

void Logger::enableConsole(bool enable) {
    std::lock_guard<std::mutex> lock(mutex_);
    console_enabled_ = enable;
}

void Logger::debug(const std::string& message) {
    log(Level::DEBUG, message);
}

void Logger::info(const std::string& message) {
    log(Level::INFO, message);
}

void Logger::warning(const std::string& message) {
    log(Level::WARNING, message);
}

void Logger::error(const std::string& message) {
    log(Level::ERROR, message);
}

void Logger::log(Level level, const std::string& message) {
    if (!shouldLog(level)) {
        return;
    }

    std::lock_guard<std::mutex> lock(mutex_);

    // Get timestamp
    auto now = std::chrono::system_clock::now();
    auto time_t_now = std::chrono::system_clock::to_time_t(now);
    auto ms = std::chrono::duration_cast<std::chrono::milliseconds>(
        now.time_since_epoch()) % 1000;

    // Get thread ID
    std::thread::id thread_id = std::this_thread::get_id();
    std::hash<std::thread::id> hasher;
    size_t thread_hash = hasher(thread_id);
    int thread_num = thread_hash % 1000;  // Simple thread number

    std::stringstream ss;
    ss << std::put_time(std::localtime(&time_t_now), "%Y-%m-%d %H:%M:%S");
    ss << '.' << std::setfill('0') << std::setw(3) << ms.count();
    ss << " [THREAD " << std::setw(3) << thread_num << "]";
    ss << " [" << levelToString(level) << "] " << message;

    std::string log_line = ss.str();

    // Output to console (thread-safe with mutex)
    if (console_enabled_) {
        if (level == Level::ERROR) {
            std::cerr << log_line << std::endl;
        } else {
            std::cout << log_line << std::endl;
        }
    }

    // Output to file (thread-safe with mutex)
    if (log_file_.is_open()) {
        log_file_ << log_line << std::endl;
        log_file_.flush();
    }
}

std::string Logger::levelToString(Level level) const {
    switch (level) {
        case Level::DEBUG:   return "DEBUG";
        case Level::INFO:    return "INFO";
        case Level::WARNING: return "WARNING";
        case Level::ERROR:   return "ERROR";
        default:             return "UNKNOWN";
    }
}

bool Logger::shouldLog(Level level) const {
    return level >= min_level_;
}
