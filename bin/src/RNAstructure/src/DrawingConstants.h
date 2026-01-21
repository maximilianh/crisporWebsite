#ifndef DRAWCONSTANTS_H
#define DRAWCONSTANTS_H

#include <sstream>
using namespace std;

/******************************************/
/* Text sizing.                           */
/******************************************/

const int GLYPHSIZE = 15;
const int TEXTSIZE = 24;
const string TEXTSIZE_STRING = "24";
const int TEXTSIZE_LEGEND = 16;
const string TEXTSIZE_LEGEND_STRING = "16";
const int TEXTSIZE_BIGGER = 26;
const string TEXTSIZE_BIGGER_STRING = "26";
const int TEXTSIZE_LEGEND_BIGGER = 18;
const string TEXTSIZE_LEGEND_BIGGER_STRING = "18";

/******************************************/
/* Miscellaneous layout properties.       */
/******************************************/

const int BORDER = 36;
const string BORDER_STRING = "36";
const int CIRCLE_RADIUS = 12;
const string CIRCLE_RADIUS_STRING = "12";
const double NUC_LABEL_HEIGHT = 24;
const string NUC_LABEL_HEIGHT_STRING = "24";
const int NUC_POS_ADJUSTMENT = 4;
const double PI = 3.14159;

/******************************************/
/* Color identification strings.          */
/******************************************/

// Generic color strings.
const string BLACK = "Black";
const string WHITE = "White";
const string GRAY = "Gray";
const string RED = "Red";
const string PINK = "Pink";
const string ORANGE = "Orange";
const string YELLOW = "Yellow";
const string LIGHT_GREEN = "Light Green";
const string GREEN = "Green";
const string LIGHT_BLUE = "Light Blue";
const string BLUE = "Blue";
const string PURPLE = "Purple";

// Inline function to get a particular type of color string.
inline string getColorString( string type, bool isSVG ) {
	if( type == WHITE ) {
		if( !isSVG ) { return "1.00 1.00 1.00"; }
		else { return "\"rgb(255,255,255)\""; }
	} else if( type == GRAY ) {
		if( !isSVG ) { return "0.67 0.67 0.67"; }
		else { return "\"rgb(171,171,171)\""; }
	} else if( type == RED ) {
		if( !isSVG ) { return "1.00 0.00 0.00"; }
		else { return "\"rgb(255,0,0)\""; }
	} else if( type == PINK ) {
		if( !isSVG ) { return "1.00 0.50 1.00"; }
		else { return "\"rgb(255,128,255)\""; }
	} else if( type == ORANGE ) {
		if( !isSVG ) { return "1.00 0.50 0.00"; }
		else { return "\"rgb(255,171,0)\""; }
	} else if( type == YELLOW ) {
		if( !isSVG ) { return "0.83 0.83 0.17"; }
		else { return "\"rgb(212,212,44)\""; }
	} else if( type == LIGHT_GREEN ) {
		if( !isSVG ) { return "0.00 1.00 0.00"; }
		else { return "\"rgb(0,255,0)\""; }
	} else if( type == GREEN ) {
		if( !isSVG ) { return "0.00 0.50 0.00"; }
		else { return "\"rgb(0,128,0)\""; }
	} else if( type == LIGHT_BLUE ) {
		if( !isSVG ) { return "0.00 0.67 1.00"; }
		else { return "\"rgb(0,171,255)\""; }
	} else if( type == BLUE ) {
		if( !isSVG ) { return "0.00 0.00 1.00"; }
		else { return "\"rgb(0,0,255)\""; }
	} else if( type == PURPLE ) {
		if( !isSVG ) { return "0.50 0.00 0.50"; }
		else { return "\"rgb(128,0,128)\""; }
	} else {
		if( !isSVG ) { return "0.00 0.00 0.00"; }
		else { return "\"rgb(0,0,0)\""; }
	}
}

/******************************************/
/* Image size constants.                  */
/******************************************/

// Maximum Postscript description length and image bounds.
const unsigned int DESC_PS = 39;
const int XBOUND_PS = 612;
const int YBOUND_PS = 792;

// Maximum SVG description length and image bounds.
const unsigned int DESC_SVG = 45;
const int XBOUND_SVG = 790;
const int YBOUND_SVG = 905;

/******************************************/
/* Drawing object template pieces.        */
/******************************************/

// Definitions of standardized variables used in structure element templates.
const string BACKGROUND = "BACKGROUND";
const string COLOR = "COLOR";
const string CONTROLX = "CONTROLX";
const string CONTROLY = "CONTROLY";
const string CURVEWEIGHT = "CURVEWEIGHT";
const string ENDX = "ENDX";
const string ENDY = "ENDY";
const string HEIGHT = "HEIGHT";
const string LINEWEIGHT = "LINEWEIGHT";
const string LOCX = "LOCX";
const string LOCY = "LOCY";
const string OUTLINE = "OUTLINE";
const string RADIUS = "RADIUS";
const string SCALEFACTOR = "SCALEFACTOR";
const string STARTX = "STARTX";
const string STARTY = "STARTY";
const string TEXTSTRING = "TEXTSTRING";
const string WIDTH = "WIDTH";
const string X1 = "X1";
const string X2 = "X2";
const string Y1 = "Y1";
const string Y2 = "Y2";

// Postscript color template.
const string COLOR_TEMPLATE_PS = "RED GREEN BLUE";

// Postscript scaling region markers. 
const string SCALE_OPEN_PS = "gsave " + SCALEFACTOR + " " + SCALEFACTOR + " scale";
const string SCALE_CLOSE_PS = "grestore";

// Postscript legend resizing markers.
const string LEGEND_RESIZE_START_PS = "[" + TEXTSIZE_LEGEND_STRING + " 0 0 -" + TEXTSIZE_LEGEND_STRING + " 0 0] /Courier-Bold sfm";
const string LEGEND_RESIZE_END_PS = "";

// Postscript syntax strings to draw specific structural elements.
const string CIRCLE_PS = LINEWEIGHT + " setlinewidth newpath " + OUTLINE + " setrgbcolor " + LOCX + " " + LOCY + " " + RADIUS + " 0 360 arc closepath gsave " + BACKGROUND + " setrgbcolor fill grestore stroke";
const string CURVE_PS = COLOR + " setrgbcolor " + CURVEWEIGHT + " setlinewidth " + X1 + " " + Y1 + " moveto " + X1 + " " + Y1 + " " + CONTROLX + " " + CONTROLY + " " + X2 + " " + Y2 + " curveto stroke";
const string LINE_PS = COLOR + " setrgbcolor " + LINEWEIGHT + " setlinewidth newpath " + STARTX + " " + STARTY + " moveto " + ENDX + " " + ENDY + " lineto closepath stroke";
const string RECTANGLE_PS = COLOR + " setrgbcolor newpath " + LOCX + " " + LOCY + " moveto 0 " + HEIGHT + " rlineto " + WIDTH + " 0 rlineto 0 -" + HEIGHT + " rlineto closepath fill";
const string TEXT_PS = LOCX + " " + LOCY + " moveto " + COLOR + " setrgbcolor (" + TEXTSTRING + ") show";

// SVG color template.
const string COLOR_TEMPLATE_SVG = "\"rgb(RED,GREEN,BLUE)\"";

// SVG scaling region markers.
const string SCALE_OPEN_SVG = "<g transform=\"scale(" + SCALEFACTOR + ")\">";
const string SCALE_CLOSE_SVG = "</g>";

// SVG legend resizing markers.
const string LEGEND_RESIZE_START_SVG = "<g font-size=\"" + TEXTSIZE_LEGEND_STRING + "\">";
const string LEGEND_RESIZE_END_SVG = "</g>";

// SVG syntax strings to draw specific structural elements.
const string CIRCLE_SVG = "<circle style=\"stroke-width:" + LINEWEIGHT + "\" cx=\"" + LOCX + "\" cy=\"" + LOCY + "\" r=\"" + RADIUS + "\" fill=" + BACKGROUND + " stroke=" + OUTLINE + "/>";
const string CURVE_SVG = "<path style=\"fill:none;stroke-width:" + CURVEWEIGHT + "\" stroke=" + COLOR + " d=\"M" + X1 + "," + Y1 + " Q" + CONTROLX + "," + CONTROLY + " " + X2 + "," + Y2 + "\"/>";
const string LINE_SVG = "<line style=\"fill:none;stroke-width:" + LINEWEIGHT + "\" stroke=" + COLOR + " x1=\"" + STARTX + "\" y1=\"" + STARTY + "\" x2=\"" + ENDX + "\" y2=\"" + ENDY + "\"/>";
const string RECTANGLE_SVG = "<rect x=\"" + LOCX + "\" y=\"" + LOCY + "\" width=\"" + WIDTH + "\" height=\"" + HEIGHT + "\" fill=" + COLOR + " stroke=" + COLOR + "/>";
const string TEXT_SVG = "<text x=\"" + LOCX + "\" y=\"" + LOCY + "\" fill=" + COLOR + " stroke=" + COLOR + ">" + TEXTSTRING + "</text>";

/******************************************/
/* File boundary marker creation.         */
/******************************************/

// Postscript file markers.
inline string createStartPS(int pagenumber=1, int pages=1) {
	stringstream startMarkerStream( stringstream::in | stringstream::out );
	startMarkerStream <<"%!PS-Adobe-3.0"<<endl
	                  << endl
	                  <<"%%Pages: "<<pages<<endl
	                  <<"%%Page: "<<pagenumber<<" "<<pagenumber<<endl
	                  << "0 " << YBOUND_PS << " translate 1 -1 scale" << endl
	                  << "/sfm { findfont exch makefont setfont } bind def" << endl
	                  << "[" << TEXTSIZE << " 0 0 " << -TEXTSIZE << " 0 0] /Courier-Bold sfm";
	return startMarkerStream.str();
}
const string END_MARKER_PS = "showpage";

// SVG file markers.
inline string createStartSVG() {
	stringstream startMarkerStream( stringstream::in | stringstream::out );
	startMarkerStream << "<?xml version=\"1.0\" encoding=\"UTF-8\" standalone=\"no\"?>" << endl
	                  << "<!DOCTYPE svg PUBLIC \"-//W3C//DTD SVG 1.1//EN\" "
	                  << "\"http://www.w3.org/Graphics/SVG/1.1/DTD/svg11.dtd\">" << endl
	                  << "<svg xmlns=\"http://www.w3.org/2000/svg\" "
	                  << "xmlns:xlink=\"http://www.w3.org/1999/xlink\" "
	                  << "xml:space=\"preserve\" font-family=\"monospace\" font-size=\"" << TEXTSIZE << "\" "
	                  << "fill=" << getColorString( WHITE, true ) << " stroke=" << getColorString( BLACK, true ) << " "
	                  << "viewBox=\"0 0 " << XBOUND_SVG << " " << YBOUND_SVG << "\">";
	return startMarkerStream.str();
}
const string END_MARKER_SVG = "</svg>";

// Options to specify how descriptions should be output in the legend of an image, ps, or other plot.
enum LegendDescriptionOutputType {
	// Use the default/existing description, which varies based on program. 
	// (For ProbabilityPlot, the default description is the input/source file name. For StructureImageHandler, the description is obtained from the structure comment. )
	DESC_USE_DEFAULT,
	// Do not include the description in the legend.
	DESC_USE_NONE,
	// Specifies that if (and only if) the description represents a path to a file (e.g. the input file),
	// the description that is shown in the legend should be only the base name of the file (i.e. without directory or extension).
	DESC_USE_FILENAME,
	// Indicates that a custom description should be used instead of the default. 
	// This assumes the custom text is provided by the caller (using a function parameter) or by a user (with a command-line argument) etc.
	DESC_USE_CUSTOM
};

struct LegendDescriptionSettings {
	LegendDescriptionOutputType outputType;
	vector<string> customDescriptions;
	LegendDescriptionSettings() { outputType=DESC_USE_DEFAULT; }
	string getDescription(const string defaultDescription, bool isFilePath=false, const int plotIndex=0) {
		switch(outputType) {
			case DESC_USE_NONE: return std::string();
			case DESC_USE_CUSTOM: 
				if (plotIndex >= 0 && plotIndex < customDescriptions.size())
					return customDescriptions[plotIndex];
				return defaultDescription;
			case DESC_USE_FILENAME:
				// if the default description represents a file, return only the file name.
				if (isFilePath)
					return getFileBaseName(defaultDescription);
				return defaultDescription;
			case DESC_USE_DEFAULT: 
			default: return defaultDescription;
		}
	}
	void parse(string optionText) {
		customDescriptions.clear();
		if (optionText == "~~"|| optionText == "~default")
			outputType=DESC_USE_DEFAULT;
		else if (optionText.empty() || optionText=="~none")
			outputType=DESC_USE_NONE;
		else if (optionText == "~file")
			outputType=DESC_USE_FILENAME;
		else if (optionText.compare(0, 5, "~list")==0) {
			outputType=DESC_USE_CUSTOM;
			stringSplit(customDescriptions, optionText.substr(6, string::npos), optionText[5]);
		} else {
			outputType=DESC_USE_CUSTOM;
			customDescriptions.push_back(optionText);
		}
	}
	static string getFileBaseName(const string& file) {
		const char slash[] = { '/', '\\' };
		int start = file.find_last_of(slash); //get last index of slash ( / or \ ), or string:npos if neither is found.
		if (start == string::npos) start = 0; else start++; // if a slash was found, start with the character AFTER the slash. otherwise start at 0.
		if (start == file.size()) return string(); //return empty because the slash was the last character.
		int end = file.rfind('.'); // get the last position of '.' (or string::npos if it is not found)
		if (end < start || end == string::npos) end = file.size(); // set the end to the full string length if the dot is before the slash. (e.g. /folder/name/with.slash/filename)
		return file.substr(start, end-start);
	}
	static string escapeBackSlashes(const string& s) {
		int ct = 0, sz=s.size();
		for(int i = 0; i < sz; i++)
			if (s[i]=='\\') ct++;
		if (ct==0) return s;
		string snew = string(sz+ct, 0);
		for(int i = 0, ct=0; i < sz; i++) {
			snew[ct++]=s[i];
			if (s[i]=='\\')
				snew[ct++]='\\'; //insert an extra backslash if the current character is a backslash.
		}
		return snew;
	}
	// Attempts to guess whether the string represents a file, for the purpose of deciding whether a figure description
	// represents a file path. There's not a rigorous way to do this, because the file need not exist.
	// But a good guess is that it is a file if it has an extension, because most of the input files to the draw/plot programs 
	// have extensions (.ct, .bracket etc). We'll define a file extension as a dot followed by 2 to 7 ASCII characters,
	// one of which MUST be a letter and none of which can be a space. 
	// This excludes things like "A Sentence."  "J.S. Name"  "A number 3.43"  "Elipses..."
	static bool guessTextIsFile(const string& s) {
		int len = s.length();
		int pos = s.rfind('.');  
		//if (pos == string::npos || pos == len-1) return false;
		//len = len - pos - 1; //the length of the extension
		//string ext = s.substr(pos+1, len); //get just the extension (not including the dot)
		int extLen = len - pos - 1;
		if (extLen < 2 || extLen > 7) return false;
		bool hasLetter = false;
		for (int i=pos+1; i < len; i++) {
			if ((s[i] >= 'a' && s[i] <= 'z') || (s[i] >= 'A' && s[i] <= 'Z')) 
				hasLetter = true;
			else if (s[i] == ' ') 
				return false;
		}
		return hasLetter;
		//string exts[] = { "ct", "ps", "svg", "fasta", "lis", "conf", "seq", "dot", "shape", "txt", "con", "pp", "dp", "bracket"};
		//const int numExts = sizeof(exts) / sizeof(string);
		//for (int i = 0; i < numExts; i++)
		//	if (ext==exts[i]) return true;
		//return false;
	}

private:
	static inline int spos(int pos, int alt) { if (pos==string::npos) return alt; return pos; }
	static void stringSplit(vector<string> &tokens, const string& text, char sep) {
		int start = 0, end = 0;
		while ((end = text.find(sep, start)) != string::npos) {
			tokens.push_back(text.substr(start, end - start));
			start = end + 1;
		}
		tokens.push_back(text.substr(start));
	}
};


#define LOGP_EPSILON 1.0E-10 // Represents a VERY high probability, close to 1. E.g. -log10(0.999999999) = 4.34E-10
#define LOGP_ZERO(val) ((val<0?-val:val)<LOGP_EPSILON?0:val)

#endif /* DRAWINGCONSTANTS_H */
