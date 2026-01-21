//! Determine the line number of the current read-position of an 
//! input stream. This process clears the current state of the 
//! stream, so it should only be used when the stream is no longer needed,
//! for example after an error occurs to show an informative error message 
//! to the user.
//! This method attempts to store the current position of an input stream,
//! then reset and rescan it to count the number of newline characters 
//! encountered up to the previous (stored) position.
//!
//! The line number returned is 1-based. 
//! A return value of 0 indicates that the stream was not opened.
unsigned long int getCurrentLine(std::istream& is);

//! A class that implments a filtering stream buffer that tracks the current
//! line number in a stream (by counting linefeeds).
//! Useful for showing informative error messages when reading files.
class LineCounterStreambuf : public std::streambuf {
    std::streambuf* source;
    std::istream* owner;
    bool atStartOfLine;
    unsigned int line;
    char buf;
protected:
    int underflow();
public:
    LineCounterStreambuf( std::streambuf* source );
    LineCounterStreambuf( std::istream& owner );
    ~LineCounterStreambuf();
    int lineNumber() const;
};
