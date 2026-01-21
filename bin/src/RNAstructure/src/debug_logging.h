/*
Implementation of the logging framework described in http://www.drdobbs.com/cpp/logging-in-c/201804215


*/

#ifndef DEBUG_LOGGING

#include <sstream>
#include <string>
#include <stdio.h>

#define DEBUG_LOGGING

/*************************************************************************
  Preprocessor Macros to control the logging behavior
**************************************************************************/
// Debug logs are disabled by default.  To enable, ENABLE_DEBUG_LOGS must
// be defined before debug_logging.h is included as a header file.  If 
// ENABLE_DEBUG_LOGS is not defined, then the macros are replaced by empty
// strings, disabling all logging functionality.
#ifdef ENABLE_DEBUG_LOGS
	// The general call to the logger.  level is one of the debugging levels
	// enumerated in LogLevel below.  After next arguments are passed to the
	// stringstream.  An example usage is:
	//
	//    DEBUG_LOG(DEBUG4, "cumulative: " << cumulative)
	//
    #define DEBUG_LOG(levelName, ...) \
    if (DebugLog::FromString(#levelName) > DebugLog::ReportingLevel() || !OutputHandler::Stream()); \
    else DebugLog().Write(DebugLog::FromString(#levelName)) << __FILE__ << "\t" << __FUNCTION__ << "\t" << __LINE__ << "\t" << __VA_ARGS__ 

	// The logger outputs to std::stderr by default.  This behavior can be changed 
	// with the following macro, which causes the logger to write to a file. Example
	// usage:
	//
	//     USE_LOG_FILE("foo.log")
	//
    #define LOG_TO_FILE(log_file) FILE* pFile = fopen(log_file, "a"); OutputHandler::Stream() = pFile

    // This macro switches back to stderr
    #define LOG_TO_STDERR() FILE* pStream = stderr; OutputHandler::Stream() = pStream
	
    // This macro changes the level of debugging messages that are output.  
	// Acceptable levelName values are::
	//   "ERROR", 
    //   "WARNING", 
    //   "INFO", 
    //   "DEBUG", 
    //   "DEBUG1", 
    //   "DEBUG2", 
    //   "DEBUG3", 
    //   "DEBUG4", 
    //   "TRACE"
	//
    #define SET_DEBUG_LEVEL(levelName) DebugLog::ReportingLevel() = DebugLog::FromString(#levelName)
    // if (logDEBUG > DebugLog::ReportingLevel() || !OutputHandler::Stream()); \
    // else DebugLog().Write(logDEBUG) << __FILE__ << "\t" << __FUNCTION__ << "\t" << __LINE__ << "\t" << "Setting log level to " << #levelName;    

    #define LOG_ERROR( ...) \
    if (logERROR > DebugLog::ReportingLevel() || !OutputHandler::Stream()); \
    else DebugLog().Write(logERROR) << __FILE__ << "\t" << __FUNCTION__ << "\t" << __LINE__ << "\t" << __VA_ARGS__ 

    #define LOG_WARN( ...) \
    if (logWARNING > DebugLog::ReportingLevel() || !OutputHandler::Stream()); \
    else DebugLog().Write(logWARNING) << __FILE__ << "\t" << __FUNCTION__ << "\t" << __LINE__ << "\t" << __VA_ARGS__ 

    #define LOG_INFO( ...) \
    if (logINFO > DebugLog::ReportingLevel() || !OutputHandler::Stream()); \
    else DebugLog().Write(logINFO) << __FILE__ << "\t" << __FUNCTION__ << "\t" << __LINE__ << "\t" << __VA_ARGS__ 

    #define LOG_DEBUG( ...) \
    if (logDEBUG > DebugLog::ReportingLevel() || !OutputHandler::Stream()); \
    else DebugLog().Write(logDEBUG) << __FILE__ << "\t" << __FUNCTION__ << "\t" << __LINE__ << "\t" << __VA_ARGS__ 

    #define LOG_DEBUG2( ...) \
    if (logDEBUG2 > DebugLog::ReportingLevel() || !OutputHandler::Stream()); \
    else DebugLog().Write(logDEBUG2) << __FILE__ << "\t" << __FUNCTION__ << "\t" << __LINE__ << "\t" << __VA_ARGS__ 

    #define LOG_DEBUG3( ...) \
    if (logDEBUG3 > DebugLog::ReportingLevel() || !OutputHandler::Stream()); \
    else DebugLog().Write(logDEBUG3) << __FILE__ << "\t" << __FUNCTION__ << "\t" << __LINE__ << "\t" << __VA_ARGS__ 

    #define LOG_DEBUG4( ...) \
    if (logDEBUG4 > DebugLog::ReportingLevel() || !OutputHandler::Stream()); \
    else DebugLog().Write(logDEBUG4) << __FILE__ << "\t" << __FUNCTION__ << "\t" << __LINE__ << "\t" << __VA_ARGS__ 

    #define LOG_TRACE( ...) \
    if (logTRACE > DebugLog::ReportingLevel() || !OutputHandler::Stream()); \
    else DebugLog().Write(logTRACE) << __FILE__ << "\t" << __FUNCTION__ << "\t" << __LINE__ << "\t" << __VA_ARGS__ 

#else

    // Disables all logging functions
    #define DEBUG_LOG(...)
    #define USE_LOG_FILE(log_file)
    #define SET_DEBUG_LEVEL(levelName)
    #define LOG_ERROR( ...)  DebugLog().Write(logERROR) << __VA_ARGS__ 
    #define LOG_WARN( ...) DebugLog().Write(logWARNING) << __VA_ARGS__
    #define LOG_INFO( ...)
    #define LOG_DEBUG( ...)
    #define LOG_DEBUG2( ...)
    #define LOG_DEBUG3( ...)
    #define LOG_DEBUG4( ...)
    #define LOG_TRACE( ...)

#endif //ENABLE_DEBUG_LOGS


/*************************************************************************
  Logging class definitions
**************************************************************************/

/*
* Supported Log Levels.
*/
enum LogLevel {
    logERROR = 0, 
    logWARNING, 
    logINFO, 
    logDEBUG, 
    logDEBUG1, 
    logDEBUG2, 
    logDEBUG3, 
    logDEBUG4, 
    logTRACE
};

/*
* Logger class.
* Inherited by DebugLog. DebugLog should be used instead of Logger
*/
template <typename T> class Logger
{
    public:
        Logger();
        virtual ~Logger();
        std::ostringstream& Write(LogLevel level = logINFO);
        static LogLevel&    ReportingLevel();
        static std::string  ToString(LogLevel level);
        static LogLevel     FromString(const std::string& level);

    protected:
        std::ostringstream  os;

    private:
        Logger(const Logger&);
        Logger& operator =(const Logger&);
};

template <typename T> Logger<T>::Logger()
{
}

// Gets the stringstream, puts the debug message level into it, and then returns the stream
template <typename T> std::ostringstream& Logger<T>::Write(LogLevel level)
{
    os << ToString(level) << "\t";
    return os;
}

/*
* The message accumated in the string stream is output when the logger instance is destroyed (deconstructed?).
*/
template <typename T> Logger<T>::~Logger()
{
    os << std::endl;
    T::Output(os.str());
}

template <typename T> LogLevel& Logger<T>::ReportingLevel()
{
    static LogLevel reportingLevel = logINFO;
    return reportingLevel;
}

/*
* Returns the Log Level name from the Log Level number.
*/
template <typename T> std::string Logger<T>::ToString(LogLevel level)
{
    static const char* const levels[] = { 
        "ERROR", 
        "WARNING", 
        "INFO", 
        "DEBUG", 
        "DEBUG1", 
        "DEBUG2", 
        "DEBUG3", 
        "DEBUG4", 
        "TRACE"
    };
    return levels[level];
}

/*
* Returns the Log Level number from the Log Level name.
*/
template <typename T> LogLevel Logger<T>::FromString(const std::string& level)
{
    if (level == "TRACE")
        return logTRACE;
    if (level == "DEBUG4")
        return logDEBUG4;
    if (level == "DEBUG3")
        return logDEBUG3;
    if (level == "DEBUG2")
        return logDEBUG2;
    if (level == "DEBUG1")
        return logDEBUG1;
    if (level == "DEBUG")
        return logDEBUG;
    if (level == "INFO")
        return logINFO;
    if (level == "WARNING")
        return logWARNING;
    if (level == "ERROR")
        return logERROR;
    Logger<T>().Write(logWARNING) << "Unknown logging level '" << level << "'. Using INFO level as default.";
    return logINFO;
}


class OutputHandler
{
    public:
        static FILE*& Stream();
        static void Output(const std::string& msg);
};

inline FILE*& OutputHandler::Stream()
{
    static FILE* pStream = stderr;
    return pStream;
}

inline void OutputHandler::Output(const std::string& msg)
{   
    FILE* pStream = Stream();
    if (!pStream) return;

    fprintf(pStream, "%s", msg.c_str());
    fflush(pStream);
}

class DebugLog : public Logger<OutputHandler> {};

#endif //DEBUG_LOGGING