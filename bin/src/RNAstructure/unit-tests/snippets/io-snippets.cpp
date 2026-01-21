// Determine the line number of the current read-position of an 
// input stream. This process clears the current state of the 
// stream, so it should only be used when the stream is no longer needed,
// for example after an error occurs to show an informative error message 
// to the user.
// This method attempts to store the current position of an input stream,
// then reset and rescan it to count the number of newline characters 
// encountered up to the previous (stored) position.
//
// The line number returned is 1-based. 
// A return value of 0 indicates that the stream was not opened.
unsigned long int getCurrentLine(std::istream& is) {
    unsigned long int line = 1;
    is.clear(); // need to clear error bits otherwise tellg returns -1.
    std::streampos originalPos = is.tellg();
    if (originalPos < 0) return 0;
    is.seekg(0);
    char c;
	// TODO: Update to using sentry and rdbuf to improve efficiency.
	// TODO: Check for old-style Mac line endings ('\r')
    while ((is.tellg() < originalPos) && is.get(c))
        if (c == '\n') ++line;
    
    // ok but if we are AT a newline, subtract one.
    if (c == '\n') --line;

    return line;
}



// A class that implments a filtering stream buffer that tracks the 
// current line number in a stream (by counting linefeeds).
// Useful for showing informative error messages when reading files.
int LineCounterStreambuf::underflow() {
	int ch = source->sbumpc();
	if (ch != EOF) {
		buf = ch;
		setg(&buf, &buf, &buf+1);
		if (atStartOfLine)
			++line;
		atStartOfLine = buf == '\n';
	}
	return ch;
}
LineCounterStreambuf::LineCounterStreambuf( std::streambuf* source )
        : source(source), owner(nullptr), atStartOfLine(true), line(0) {}
LineCounterStreambuf::LineCounterStreambuf( std::istream& owner )
        : source(owner.rdbuf()), owner(&owner), atStartOfLine(true), line(0) {
	owner.rdbuf(this);
}
LineCounterStreambuf::~LineCounterStreambuf() {
	if (owner!=NULL) owner->rdbuf(source);
}
int LineCounterStreambuf::lineNumber() const {
	return line;
}