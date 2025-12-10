#pragma once

#include <string>
#include <iostream>
#include <fstream>
#include <mutex>
#include <memory>

class Logger {
public:
    enum class Level {
        DEBUG,
        INFO,
        WARNING,
        ERROR
    };

    static Logger& getInstance();

    void setLogFile(const std::string& filename);
    void setLevel(Level level);
    void enableConsole(bool enable);

    void debug(const std::string& message);
    void info(const std::string& message);
    void warning(const std::string& message);
    void error(const std::string& message);

    // Template for formatting
    template<typename... Args>
    void info(const std::string& format, Args... args) {
        log(Level::INFO, format, args...);
    }

    template<typename... Args>
    void debug(const std::string& format, Args... args) {
        log(Level::DEBUG, format, args...);
    }

    template<typename... Args>
    void warning(const std::string& format, Args... args) {
        log(Level::WARNING, format, args...);
    }

    template<typename... Args>
    void error(const std::string& format, Args... args) {
        log(Level::ERROR, format, args...);
    }

private:
    Logger();
    ~Logger();
    Logger(const Logger&) = delete;
    Logger& operator=(const Logger&) = delete;

    void log(Level level, const std::string& message);

    template<typename... Args>
    void log(Level level, const std::string& format, Args... args) {
        char buffer[1024];
        snprintf(buffer, sizeof(buffer), format.c_str(), args...);
        log(level, std::string(buffer));
    }

    std::string levelToString(Level level) const;
    bool shouldLog(Level level) const;

    Level min_level_;
    bool console_enabled_;
    std::ofstream log_file_;
    std::mutex mutex_;
};

// Convenience macros
#define LOG_DEBUG(msg) Logger::getInstance().debug(msg)
#define LOG_INFO(msg) Logger::getInstance().info(msg)
#define LOG_WARNING(msg) Logger::getInstance().warning(msg)
#define LOG_ERROR(msg) Logger::getInstance().error(msg)
