/*
 * An implementation for a class that takes two CT files, then compares their pairings and outputs the data graphically.
 *
 * (c) 2010 Mathews Lab, University of Rochester Medical Center.
 * Written by Jessica S. Reuter.
 */

#include "StructureComparedImageHandler.h"

///////////////////////////////////////////////////////////////////////////////
// Constructor.
///////////////////////////////////////////////////////////////////////////////
StructureComparedImageHandler::StructureComparedImageHandler() {

	numPredictedStructures = 0;
}

///////////////////////////////////////////////////////////////////////////////
// Write additional information to show comparison statistics.
///////////////////////////////////////////////////////////////////////////////
void StructureComparedImageHandler::addComparisonData( string predicted, string accepted, int number, bool exact, bool isSVG, bool filenames ) {

	// Initialize an index used to handle description sanitization and conversion template indices later.
	size_t index = string::npos;

	// Get the predicted and accepted structures. Error checking isn't necessary because it's been done before.
	// Also get the length of each strand, since they're the same.
	RNA* predictedStrand = new RNA( predicted.c_str(), FILE_CT );
	RNA* acceptedStrand = new RNA( accepted.c_str(), FILE_CT );
	int length = predictedStrand->GetSequenceLength();

	// Calculate sensitivity.
	int pairs1 = 0;
	int score1 = 0;
	stringstream sensitivity( stringstream::in | stringstream::out );
	scorer( acceptedStrand->GetStructure(), predictedStrand->GetStructure(), &score1, &pairs1, number, exact );
	double percent1 = ( ( (double)score1 ) / ( (double)pairs1 ) ) * 100;
	sensitivity << "Sensitivity: " << score1 << " / " << pairs1 << " = " << fixed << setprecision( 2 ) << percent1 << "%";

	// Calculate PPV.
	int pairs2 = 0;
	int score2 = 0;
	stringstream ppv( stringstream::in | stringstream::out );
	scorerppv( acceptedStrand->GetStructure(), predictedStrand->GetStructure(), &score2, &pairs2, number, exact );
	double percent2 = ( ( (double)score2 ) / ( (double)pairs2 ) ) * 100;
	ppv << "PPV: " << score2 << " / " << pairs2 << " = " << fixed << setprecision( 2 ) << percent2 << "%";

	// Calculate the total number of pairs and the total number of pseudoknotted pairs in the predicted structure.
	int predictedPairsNum = 0, predictedPseudoknottedPairsNum = 0, predictedAfterBreakageNum = 0;
	for( int i = 1; i <= length; i++ ) {
		if( predictedStrand->GetPair( i, number ) > i ) { predictedPairsNum++; }
	}
	predictedStrand->BreakPseudoknot( false, number );
	for( int i = 1; i <= length; i++ ) {
		if( predictedStrand->GetPair( i, number ) > i ) { predictedAfterBreakageNum++; }
	}
	predictedPseudoknottedPairsNum = predictedPairsNum - predictedAfterBreakageNum;

	// Calculate the total number of pairs and the total number of pseudoknotted pairs in the accepted structure.
	int acceptedPairsNum = 0, acceptedPseudoknottedPairsNum = 0, acceptedAfterBreakageNum = 0;
	for( int i = 1; i <= length; i++ ) {
		if( acceptedStrand->GetPair( i ) > i ) { acceptedPairsNum++; }
	}
	acceptedStrand->BreakPseudoknot( false );
	for( int i = 1; i <= length; i++ ) {
		if( acceptedStrand->GetPair( i ) > i ) { acceptedAfterBreakageNum++; }
	}
	acceptedPseudoknottedPairsNum = acceptedPairsNum - acceptedAfterBreakageNum;

	// Get the predicted structure description.
	string predictedDesc = predictedStrand->GetCommentString( number );
	if( true ) {
		const size_t whiteIndex1 = predictedDesc.find_first_not_of( " \n\r\t" );
		if( whiteIndex1 != string::npos ) { predictedDesc = predictedDesc.substr( whiteIndex1 ); }
		const size_t whiteIndex2 = predictedDesc.find_last_not_of( " \n\r\t" );
		if( whiteIndex2 != string::npos ) { predictedDesc = predictedDesc.substr( 0, whiteIndex2 + 1 ); }
		while( (index = predictedDesc.find( '(', 0 )) != string::npos ) { predictedDesc = predictedDesc.erase( index, 1 ); }
		while( (index = predictedDesc.find( ')', 0 )) != string::npos ) { predictedDesc = predictedDesc.erase( index, 1 ); }
		while( (index = predictedDesc.find( '<', 0 )) != string::npos ) { predictedDesc = predictedDesc.erase( index, 1 ); }
		while( (index = predictedDesc.find( '>', 0 )) != string::npos ) { predictedDesc = predictedDesc.erase( index, 1 ); }
	}

	// Get the accepted structure description.
	string acceptedDesc = acceptedStrand->GetCommentString();
	if( true ) {
		const size_t whiteIndex1 = acceptedDesc.find_first_not_of( " \n\r\t" );
		if( whiteIndex1 != string::npos ) { acceptedDesc = acceptedDesc.substr( whiteIndex1 ); }
		const size_t whiteIndex2 = acceptedDesc.find_last_not_of( " \n\r\t" );
		if( whiteIndex2 != string::npos ) { acceptedDesc = acceptedDesc.substr( 0, whiteIndex2 + 1 ); }
		while( (index = acceptedDesc.find( '(', 0 )) != string::npos ) { acceptedDesc = acceptedDesc.erase( index, 1 ); }
		while( (index = acceptedDesc.find( ')', 0 )) != string::npos ) { acceptedDesc = acceptedDesc.erase( index, 1 ); }
		while( (index = acceptedDesc.find( '<', 0 )) != string::npos ) { acceptedDesc = acceptedDesc.erase( index, 1 ); }
		while( (index = acceptedDesc.find( '>', 0 )) != string::npos ) { acceptedDesc = acceptedDesc.erase( index, 1 ); }
	}

	// Delete the strand data structures.
	delete predictedStrand;
	delete acceptedStrand;

	// Determine the added text size and legend border.
	int addedTextSize = ( !isSVG ) ? ADDED_TEXT_SIZE_PS : ADDED_TEXT_SIZE_SVG;
	int legendBorder = 10;

	// Determine the maximum bounds and the glyph size.
	int maxX = ( !isSVG ) ? XBOUND_PS : XBOUND_SVG;
	int maxY = ( !isSVG ) ? YBOUND_PS : YBOUND_SVG;
	double glyphSize = ( !isSVG ) ? ADDED_GLYPH_SIZE_PS : ADDED_GLYPH_SIZE_SVG;

	// Convert the border value to a string for alignment purposes.
	stringstream borderStream( stringstream::in | stringstream::out );
	borderStream << legendBorder;
	string borderString = borderStream.str();

	// Clear any previous comparison data.
	extras.clear();

	// Open the scaled text area.
	extras.push_back( ( !isSVG ) ? ADDITIONAL_TEXT_OPEN_PS : ADDITIONAL_TEXT_OPEN_SVG );

	// If a legend should be copied into the class, copy it.
	// The superclass part of the writeImageFile method isn't used here because the legend is placed in the upper left hand corner of the image,
	// rather than below it.
	if( legend.size() > 0 ) {

		// If the legend hasn't been copied yet, do so.
		if( legendCompared.size() == 0 ) {
			unsigned int legendPlusOne = legend.size() + 1;
			for( unsigned int i = 1; i <= legend.size(); i++ ) {

				// Determine the Y location for this entry string.
				double location = legendBorder + ( addedTextSize * i );
				stringstream entryStream( stringstream::in | stringstream::out );
				entryStream << location;
				string entryYString = entryStream.str();

				// Build the next legend entry string and add it to the legend vector in this class.
				string colorKey = legendColors[i-1];
				string text = ( isSVG ) ? TEXT_SVG : TEXT_PS;
				while( ( index = text.find( COLOR ) ) != string::npos ) { text = text.replace( index, COLOR.size(), getColorString( colorKey, isSVG ) ); }
				while( ( index = text.find( LOCX ) ) != string::npos ) { text = text.replace( index, LOCX.size(), borderString ); }
				while( ( index = text.find( LOCY ) ) != string::npos ) { text = text.replace( index, LOCY.size(), entryYString ); }
				while( ( index = text.find( TEXTSTRING ) ) != string::npos ) { text = text.replace( index, TEXTSTRING.size(), legend[i-1] ); }
				legendCompared.push_back( text );
			}

			// Clear the legend from the superclass.
			legend.clear();
		}
	}

	// If a copied legend should be written into the class, do so.
	if( legendCompared.size() > 0 ) {
		for( unsigned int i = 1; i <= legendCompared.size(); i++ ) { extras.push_back( legendCompared[i-1] ); }
	}

	// Calculate where the sensitivity string should go and place it in the image.
	if( true ) {

		// Get the sensitivity string.
		string sensitivityString = sensitivity.str();

		// Determine where the X coordinate of the string should be.
		double xLocation = ( maxX - BORDER ) - ( sensitivityString.length() * glyphSize );
		stringstream xStream( stringstream::in | stringstream::out );
		xStream << xLocation;

		// Determine where the Y coordinate of the string should be.
		stringstream yStream( stringstream::in | stringstream::out );
		yStream << ( legendBorder + addedTextSize );

		// Add the sensitivity string to the image.
		string text = ( isSVG ) ? TEXT_SVG : TEXT_PS;
		while( ( index = text.find( COLOR ) ) != string::npos ) { text = text.replace( index, COLOR.size(), getColorString( BLACK, isSVG ) ); }
		while( ( index = text.find( LOCX ) ) != string::npos ) { text = text.replace( index, LOCX.size(), xStream.str() ); }
		while( ( index = text.find( LOCY ) ) != string::npos ) { text = text.replace( index, LOCY.size(), yStream.str() ); }
		while( ( index = text.find( TEXTSTRING ) ) != string::npos ) { text = text.replace( index, TEXTSTRING.size(), sensitivityString ); }
		extras.push_back( text );
	}

	// Calculate where the PPV string should go and place it in the image.
	if( true ) {

		// Get the PPV string.
		string ppvString = ppv.str();

		// Determine where the X coordinate of the string should be.
		double adjust = ( !isSVG ) ? 1.5 : 1;
		double xLocation = ( maxX - BORDER ) - ( ( ppvString.length() - adjust ) * glyphSize );
		stringstream xStream( stringstream::in | stringstream::out );
		xStream << xLocation;

		// Determine where the Y coordinate of the string should be.
		stringstream yStream( stringstream::in | stringstream::out );
		yStream << ( legendBorder + ( addedTextSize * 2 ) );

		// Add the PPV string to the image.
		string text = ( isSVG ) ? TEXT_SVG : TEXT_PS;
		while( ( index = text.find( COLOR ) ) != string::npos ) { text = text.replace( index, COLOR.size(), getColorString( BLACK, isSVG ) ); }
		while( ( index = text.find( LOCX ) ) != string::npos ) { text = text.replace( index, LOCX.size(), xStream.str() ); }
		while( ( index = text.find( LOCY ) ) != string::npos ) { text = text.replace( index, LOCY.size(), yStream.str() ); }
		while( ( index = text.find( TEXTSTRING ) ) != string::npos ) { text = text.replace( index, TEXTSTRING.size(), ppvString ); }
		extras.push_back( text );
	}

	// Add the predicted and accepted structure text data strings.
	if( true ) {

		/// Determine the Y location of each piece of text, and convert those locations into strings.
		string y1, y2, y3;
		int xBound = ( !isSVG ) ? XBOUND_PS : XBOUND_SVG;
		int addedTextSize = ( !isSVG ) ? ADDED_TEXT_SIZE_PS : ADDED_TEXT_SIZE_SVG;
		stringstream locationStream( stringstream::in | stringstream::out );
		locationStream << ( xBound - (addedTextSize * 4) ) << " " << ( xBound - (addedTextSize * 3) ) << " " << ( xBound - (addedTextSize * 2) );
		locationStream >> y1;
		locationStream >> y2;
		locationStream >> y3;

		// Add the predicted structure text data strings.
		if( true ) {

			// Add the label.		
			string text1 = ( isSVG ) ? TEXT_SVG : TEXT_PS;
			while( ( index = text1.find( COLOR ) ) != string::npos ) { text1 = text1.replace( index, COLOR.size(), getColorString( BLACK, isSVG ) ); }
			while( ( index = text1.find( LOCX ) ) != string::npos ) { text1 = text1.replace( index, LOCX.size(), borderString ); }
			while( ( index = text1.find( LOCY ) ) != string::npos ) { text1 = text1.replace( index, LOCY.size(), y1 ); }
			while( ( index = text1.find( TEXTSTRING ) ) != string::npos ) { text1 = text1.replace( index, TEXTSTRING.size(), LABEL_PREDICTED ); }
			extras.push_back( text1 );

			// Add the number of pairs string.
			stringstream labelStream2( stringstream::in | stringstream::out );
			labelStream2 << LABEL_PAIRS << " " << predictedPairsNum;
			string label2 = labelStream2.str();
			string text2 = ( isSVG ) ? TEXT_SVG : TEXT_PS;
			while( ( index = text2.find( COLOR ) ) != string::npos ) { text2 = text2.replace( index, COLOR.size(), getColorString( BLACK, isSVG ) ); }
			while( ( index = text2.find( LOCX ) ) != string::npos ) { text2 = text2.replace( index, LOCX.size(), borderString ); }
			while( ( index = text2.find( LOCY ) ) != string::npos ) { text2 = text2.replace( index, LOCY.size(), y2 ); }
			while( ( index = text2.find( TEXTSTRING ) ) != string::npos ) { text2 = text2.replace( index, TEXTSTRING.size(), label2 ); }
			extras.push_back( text2 );

			// Add the number of pseudoknotted pairs string.
			stringstream labelStream3( stringstream::in | stringstream::out );
			labelStream3 << LABEL_PAIRS_PSEUDOKNOTTED << " " << predictedPseudoknottedPairsNum;
			string label3 = labelStream3.str();
			string text3 = ( isSVG ) ? TEXT_SVG : TEXT_PS;
			while( ( index = text3.find( COLOR ) ) != string::npos ) { text3 = text3.replace( index, COLOR.size(), getColorString( BLACK, isSVG ) ); }
			while( ( index = text3.find( LOCX ) ) != string::npos ) { text3 = text3.replace( index, LOCX.size(), borderString ); }
			while( ( index = text3.find( LOCY ) ) != string::npos ) { text3 = text3.replace( index, LOCY.size(), y3 ); }
			while( ( index = text3.find( TEXTSTRING ) ) != string::npos ) { text3 = text3.replace( index, TEXTSTRING.size(), label3 ); }
			extras.push_back( text3 );
		}

		// Add the accepted structure text data strings.
		if( true ) {

			// Set the adjustments to move the accepted strings right.
			int rightAdjust = 6;
			int rightAdjust2 = ( !isSVG ) ? 5 : 3;

			// Add the label.
			double xLocation1 = ( maxX - BORDER ) - ( LABEL_ACCEPTED.length() * glyphSize ) + ( 2 * rightAdjust ) + rightAdjust2;
			stringstream xStream1( stringstream::in | stringstream::out );
			xStream1 << xLocation1;		
			string text1 = ( isSVG ) ? TEXT_SVG : TEXT_PS;
			while( ( index = text1.find( COLOR ) ) != string::npos ) { text1 = text1.replace( index, COLOR.size(), getColorString( BLACK, isSVG ) ); }
			while( ( index = text1.find( LOCX ) ) != string::npos ) { text1 = text1.replace( index, LOCX.size(), xStream1.str() ); }
			while( ( index = text1.find( LOCY ) ) != string::npos ) { text1 = text1.replace( index, LOCY.size(), y1 ); }
			while( ( index = text1.find( TEXTSTRING ) ) != string::npos ) { text1 = text1.replace( index, TEXTSTRING.size(), LABEL_ACCEPTED ); }
			extras.push_back( text1 );

			// Add the number of pairs string.	
			stringstream labelStream2( stringstream::in | stringstream::out );
			labelStream2 << LABEL_PAIRS << " " << acceptedPairsNum;
			string label2 = labelStream2.str();
			double xLocation2 = ( maxX - BORDER ) - ( label2.length() * glyphSize ) + ( 2 * rightAdjust ) + rightAdjust2;
			stringstream xStream2( stringstream::in | stringstream::out );
			xStream2 << xLocation2;
			string text2 = ( isSVG ) ? TEXT_SVG : TEXT_PS;
			while( ( index = text2.find( COLOR ) ) != string::npos ) { text2 = text2.replace( index, COLOR.size(), getColorString( BLACK, isSVG ) ); }
			while( ( index = text2.find( LOCX ) ) != string::npos ) { text2 = text2.replace( index, LOCX.size(), xStream2.str() ); }
			while( ( index = text2.find( LOCY ) ) != string::npos ) { text2 = text2.replace( index, LOCY.size(), y2 ); }
			while( ( index = text2.find( TEXTSTRING ) ) != string::npos ) { text2 = text2.replace( index, TEXTSTRING.size(), label2 ); }
			extras.push_back( text2 );

			// Add the number of pseudoknotted pairs string.
			stringstream labelStream3( stringstream::in | stringstream::out );
			labelStream3 << LABEL_PAIRS_PSEUDOKNOTTED << " " << acceptedPseudoknottedPairsNum;
			string label3 = labelStream3.str();
			double xLocation3 = ( maxX - BORDER ) - ( label3.length() * glyphSize ) + rightAdjust;
			stringstream xStream3( stringstream::in | stringstream::out );
			xStream3 << xLocation3;
			string text3 = ( isSVG ) ? TEXT_SVG : TEXT_PS;
			while( ( index = text3.find( COLOR ) ) != string::npos ) { text3 = text3.replace( index, COLOR.size(), getColorString( BLACK, isSVG ) ); }
			while( ( index = text3.find( LOCX ) ) != string::npos ) { text3 = text3.replace( index, LOCX.size(), xStream3.str() ); }
			while( ( index = text3.find( LOCY ) ) != string::npos ) { text3 = text3.replace( index, LOCY.size(), y3 ); }
			while( ( index = text3.find( TEXTSTRING ) ) != string::npos ) { text3 = text3.replace( index, TEXTSTRING.size(), label3 ); }
			extras.push_back( text3 );
		}
	}

	// Add the pair colors legend.
	if( true ) {

		// Add the key entry for pairs found in both structures.
		stringstream labelStream1( stringstream::in | stringstream::out );
		labelStream1 << "Pair present in both predicted and accepted structure (" << colorScheme[0] << ").";
		string label1 = labelStream1.str();
		int move1 = ( filenames ) ? 9 : 6;
		double yLocation1 = ( maxY - legendBorder ) - ( addedTextSize * move1 );
		stringstream yStream1( stringstream::in | stringstream::out );
		yStream1 << yLocation1;
		string text1 = ( isSVG ) ? TEXT_SVG : TEXT_PS;
		while( ( index = text1.find( COLOR ) ) != string::npos ) { text1 = text1.replace( index, COLOR.size(), getColorString( colorScheme[0], isSVG ) ); }
		while( ( index = text1.find( LOCX ) ) != string::npos ) { text1 = text1.replace( index, LOCX.size(), borderString ); }
		while( ( index = text1.find( LOCY ) ) != string::npos ) { text1 = text1.replace( index, LOCY.size(), yStream1.str() ); }
		while( ( index = text1.find( TEXTSTRING ) ) != string::npos ) { text1 = text1.replace( index, TEXTSTRING.size(), label1 ); }
		extras.push_back( text1 );

		// Add the key entry for pairs found in the predicted structure only.
		stringstream labelStream2( stringstream::in | stringstream::out );
		labelStream2 << "Pair present in predicted structure only (" << colorScheme[1] << ").";
		string label2 = labelStream2.str();
		int move2 = ( filenames ) ? 8 : 5;
		double yLocation2 = ( maxY - legendBorder ) - ( addedTextSize * move2 );
		stringstream yStream2( stringstream::in | stringstream::out );
		yStream2 << yLocation2;
		string text2 = ( isSVG ) ? TEXT_SVG : TEXT_PS;
		while( ( index = text2.find( COLOR ) ) != string::npos ) { text2 = text2.replace( index, COLOR.size(), getColorString( colorScheme[1], isSVG ) ); }
		while( ( index = text2.find( LOCX ) ) != string::npos ) { text2 = text2.replace( index, LOCX.size(), borderString ); }
		while( ( index = text2.find( LOCY ) ) != string::npos ) { text2 = text2.replace( index, LOCY.size(), yStream2.str() ); }
		while( ( index = text2.find( TEXTSTRING ) ) != string::npos ) { text2 = text2.replace( index, TEXTSTRING.size(), label2 ); }
		extras.push_back( text2 );

		// Add the key entry for pairs found in the accepted structure only.
		stringstream labelStream3( stringstream::in | stringstream::out );
		labelStream3 << "Pair present in accepted structure only (" << colorScheme[2] << ").";
		string label3 = labelStream3.str();
		int move3 = ( filenames ) ? 7 : 4;
		double yLocation3 = ( maxY - legendBorder ) - ( addedTextSize * move3 );
		stringstream yStream3( stringstream::in | stringstream::out );
		yStream3 << yLocation3;
		string text3 = ( isSVG ) ? TEXT_SVG : TEXT_PS;
		while( ( index = text3.find( COLOR ) ) != string::npos ) { text3 = text3.replace( index, COLOR.size(), getColorString( colorScheme[2], isSVG ) ); }
		while( ( index = text3.find( LOCX ) ) != string::npos ) { text3 = text3.replace( index, LOCX.size(), borderString ); }
		while( ( index = text3.find( LOCY ) ) != string::npos ) { text3 = text3.replace( index, LOCY.size(), yStream3.str() ); }
		while( ( index = text3.find( TEXTSTRING ) ) != string::npos ) { text3 = text3.replace( index, TEXTSTRING.size(), label3 ); }
		extras.push_back( text3 );
	}

	// If applicable, add the file names which hold the predicted and accepted structures.
	if( filenames ) {

		// Add the file name of the predicted structure.
		stringstream labelStream1( stringstream::in | stringstream::out );
		labelStream1 << LABEL_PREDICTED_FILE << " " << predicted;
		string label1 = labelStream1.str();
		double yLocation1 = ( maxY - legendBorder ) - ( addedTextSize * 5 );
		stringstream yStream1( stringstream::in | stringstream::out );
		yStream1 << yLocation1;
		string text1 = ( isSVG ) ? TEXT_SVG : TEXT_PS;
		while( ( index = text1.find( COLOR ) ) != string::npos ) { text1 = text1.replace( index, COLOR.size(), getColorString( BLACK, isSVG ) ); }
		while( ( index = text1.find( LOCX ) ) != string::npos ) { text1 = text1.replace( index, LOCX.size(), borderString ); }
		while( ( index = text1.find( LOCY ) ) != string::npos ) { text1 = text1.replace( index, LOCY.size(), yStream1.str() ); }
		while( ( index = text1.find( TEXTSTRING ) ) != string::npos ) { text1 = text1.replace( index, TEXTSTRING.size(), label1 ); }
		extras.push_back( text1 );

		// Add the file name of the accepted structure.
		stringstream labelStream2( stringstream::in | stringstream::out );
		labelStream2 << LABEL_ACCEPTED_FILE << "  " << accepted;
		string label2 = labelStream2.str();
		double yLocation2 = ( maxY - legendBorder ) - ( addedTextSize * 4 );
		stringstream yStream2( stringstream::in | stringstream::out );
		yStream2 << yLocation2;
		string text2 = ( isSVG ) ? TEXT_SVG : TEXT_PS;
		while( ( index = text2.find( COLOR ) ) != string::npos ) { text2 = text2.replace( index, COLOR.size(), getColorString( BLACK, isSVG ) ); }
		while( ( index = text2.find( LOCX ) ) != string::npos ) { text2 = text2.replace( index, LOCX.size(), borderString ); }
		while( ( index = text2.find( LOCY ) ) != string::npos ) { text2 = text2.replace( index, LOCY.size(), yStream2.str() ); }
		while( ( index = text2.find( TEXTSTRING ) ) != string::npos ) { text2 = text2.replace( index, TEXTSTRING.size(), label2 ); }
		extras.push_back( text2 );
	}

	// Add the descriptions of the predicted and accepted structures.
	if( true ) {

		// Add the description of the predicted structure.
		stringstream labelStream1( stringstream::in | stringstream::out );
		labelStream1 << LABEL_PREDICTED_DESC << " " << predictedDesc;
		string label1 = labelStream1.str();
		double yLocation1 = ( maxY - legendBorder ) - ( addedTextSize * 2 );
		stringstream yStream1( stringstream::in | stringstream::out );
		yStream1 << yLocation1;
		string text1 = ( isSVG ) ? TEXT_SVG : TEXT_PS;
		while( ( index = text1.find( COLOR ) ) != string::npos ) { text1 = text1.replace( index, COLOR.size(), getColorString( BLACK, isSVG ) ); }
		while( ( index = text1.find( LOCX ) ) != string::npos ) { text1 = text1.replace( index, LOCX.size(), borderString ); }
		while( ( index = text1.find( LOCY ) ) != string::npos ) { text1 = text1.replace( index, LOCY.size(), yStream1.str() ); }
		while( ( index = text1.find( TEXTSTRING ) ) != string::npos ) { text1 = text1.replace( index, TEXTSTRING.size(), label1 ); }
		extras.push_back( text1 );

		// Add the description of the accepted structure.
		stringstream labelStream2( stringstream::in | stringstream::out );
		labelStream2 << LABEL_ACCEPTED_DESC << "  " << acceptedDesc;
		string label2 = labelStream2.str();
		double yLocation2 = ( maxY - legendBorder ) - addedTextSize;
		stringstream yStream2( stringstream::in | stringstream::out );
		yStream2 << yLocation2;
		string text2 = ( isSVG ) ? TEXT_SVG : TEXT_PS;
		while( ( index = text2.find( COLOR ) ) != string::npos ) { text2 = text2.replace( index, COLOR.size(), getColorString( BLACK, isSVG ) ); }
		while( ( index = text2.find( LOCX ) ) != string::npos ) { text2 = text2.replace( index, LOCX.size(), borderString ); }
		while( ( index = text2.find( LOCY ) ) != string::npos ) { text2 = text2.replace( index, LOCY.size(), yStream2.str() ); }
		while( ( index = text2.find( TEXTSTRING ) ) != string::npos ) { text2 = text2.replace( index, TEXTSTRING.size(), label2 ); }
		extras.push_back( text2 );
	}

	// Close the scaled text area.
	extras.push_back( ( !isSVG ) ? ADDITIONAL_TEXT_CLOSE_PS : ADDITIONAL_TEXT_CLOSE_SVG );
}

///////////////////////////////////////////////////////////////////////////////
// Get the number of predicted structures compared.
///////////////////////////////////////////////////////////////////////////////
int StructureComparedImageHandler::getNumPredictedStructures() {

	return numPredictedStructures;
}

///////////////////////////////////////////////////////////////////////////////
// Merge structure data into a single combined structure.
///////////////////////////////////////////////////////////////////////////////
void StructureComparedImageHandler::overlayStructures() {

	// Clear any previously overlaid pairs.
	pairs.clear();

	// Organize the pairs by which structure they're in.
	vector<string> both, predictedOnly, acceptedOnly;
	set_intersection( acceptedPairs.begin(), acceptedPairs.end(), predictedPairs.begin(), predictedPairs.end(), back_inserter( both ) );
	set_difference( predictedPairs.begin(), predictedPairs.end(), acceptedPairs.begin(), acceptedPairs.end(), back_inserter( predictedOnly ) );
	set_difference( acceptedPairs.begin(), acceptedPairs.end(), predictedPairs.begin(), predictedPairs.end(), back_inserter( acceptedOnly ) );

	// Colorize each pair of each pair group, then add each pair to the main pairs group.
	for( int i = 1; i <= acceptedOnly.size(); i++ ) { pairs.push_back( acceptedOnly[i-1] + " " + colorScheme[2] ); }
	for( int i = 1; i <= predictedOnly.size(); i++ ) { pairs.push_back( predictedOnly[i-1] + " " + colorScheme[1] ); }
	for( int i = 1; i <= both.size(); i++ ) { pairs.push_back( both[i-1] + " " + colorScheme[0] ); }
}

///////////////////////////////////////////////////////////////////////////////
// Read data for the accepted structure without any annotation.
///////////////////////////////////////////////////////////////////////////////
string StructureComparedImageHandler::readAcceptedStructure( string accepted ) {

	// Create a variable to track calculation results, and clear any already read pairs.
	string result = "";
	acceptedPairs.clear();

	// Read the accepted structure.
	result = readCircular( accepted, 1 );

	// Move the accepted pair data to the vector in this class.
	// Afterward, sort the accepted pair vector and clear the superclass pair vector.
	if( result == "" ) {
		int numPairs = pairs.size();
		for( int i = 1; i <= pairs.size(); i++ ) { acceptedPairs.push_back( pairs[i-1] ); }
		sort( acceptedPairs.begin(), acceptedPairs.end() );
		pairs.clear();
	}

	// Set the accepted description and clear the superclass description.
	if( result == "" ) {
		descriptionAccepted = description;
		description = "";
	}

	// Return the results string.
	return result;
}

///////////////////////////////////////////////////////////////////////////////
// Read data for the accepted structure with probability annotation.
///////////////////////////////////////////////////////////////////////////////
string StructureComparedImageHandler::readAcceptedStructureProbability( string accepted, string probabilityFile, bool text ) {

	// Create a variable to track calculation results, and clear any already read pairs.
	string result = "";
	acceptedPairs.clear();

	// Read the accepted structure.
	result = readCircular( accepted, 1 );

	// Read the annotation.
	if( result == "" ) { result = addAnnotationProbability( probabilityFile, text ); }

	// Move the accepted pair data to the vector in this class.
	// Afterward, sort the accepted pair vector and clear the superclass pair vector.
	if( result == "" ) {
		int numPairs = pairs.size();
		for( int i = 1; i <= pairs.size(); i++ ) { acceptedPairs.push_back( pairs[i-1] ); }
		sort( acceptedPairs.begin(), acceptedPairs.end() );
		pairs.clear();
	}

	// Set the accepted description and clear the superclass description.
	if( result == "" ) {
		descriptionAccepted = description;
		description = "";
	}

	// Return the results string.
	return result;
}

///////////////////////////////////////////////////////////////////////////////
// Read data for the predicted structure.
///////////////////////////////////////////////////////////////////////////////
string StructureComparedImageHandler::readPredictedStructure( string predicted, int number ) {

	// Create a variable to track calculation results, and clear any already read pairs.
	string result = "";
	predictedPairs.clear();

	// Read the predicted structure.
	result = readCircular( predicted, number );

	// Move the predicted pair data to the vector in this class.
	// Afterward, sort the predicted pair vector and clear the superclass pair vector.
	if( result == "" ) {
		int numPairs = pairs.size();
		for( int i = 1; i <= pairs.size(); i++ ) { predictedPairs.push_back( pairs[i-1] ); }
		sort( predictedPairs.begin(), predictedPairs.end() );
		pairs.clear();
	}

	// Set the predicted description and clear the superclass description.
	if( result == "" ) {
		descriptionPredicted = description;
		description = "";
	}

	// Return the results string.
	return result;
}

///////////////////////////////////////////////////////////////////////////////
// Read data for the predicted structure with probability annotation.
///////////////////////////////////////////////////////////////////////////////
string StructureComparedImageHandler::readPredictedStructureProbability( string predicted, int number, string probabilityFile, bool text ) {

	// Create a variable to track calculation results, and clear any already read pairs.
	string result = "";
	predictedPairs.clear();

	// Read the predicted structure.
	result = readCircular( predicted, number );

	// Read the annotation.
	if( result == "" ) { result = addAnnotationProbability( probabilityFile, text ); }

	// Move the predicted pair data to the vector in this class.
	// Afterward, sort the predicted pair vector and clear the superclass pair vector.
	if( result == "" ) {
		int numPairs = pairs.size();
		for( int i = 1; i <= pairs.size(); i++ ) { predictedPairs.push_back( pairs[i-1] ); }
		sort( predictedPairs.begin(), predictedPairs.end() );
		pairs.clear();
	}

	// Set the predicted description and clear the superclass description.
	if( result == "" ) {
		descriptionPredicted = description;
		description = "";
	}

	// Return the results string.
	return result;
}

///////////////////////////////////////////////////////////////////////////////
// Read data for the predicted structure with SHAPE annotation.
///////////////////////////////////////////////////////////////////////////////
string StructureComparedImageHandler::readPredictedStructureSHAPE( string predicted, int number, string SHAPEFile ) {

	// Create a variable to track calculation results, and clear any already read pairs.
	string result = "";
	predictedPairs.clear();

	// Read the predicted structure.
	result = readCircular( predicted, number );

	// Read the annotation, if it hasn't been read already.
	if( result == "" ) { result = addAnnotationSHAPE( SHAPEFile ); }

	// Move the predicted pair data to the vector in this class.
	// Afterward, sort the predicted pair vector and clear the superclass pair vector.
	if( result == "" ) {
		int numPairs = pairs.size();
		for( int i = 1; i <= pairs.size(); i++ ) { predictedPairs.push_back( pairs[i-1] ); }
		sort( predictedPairs.begin(), predictedPairs.end() );
		pairs.clear();
	}

	// Set the predicted description and clear the superclass description.
	if( result == "" ) {
		descriptionPredicted = description;
		description = "";
	}

	// Return the results string.
	return result;
}

///////////////////////////////////////////////////////////////////////////////
// Set the color scheme for the comparison.
///////////////////////////////////////////////////////////////////////////////
void StructureComparedImageHandler::setColorScheme( bool alternative ) {
	colorScheme.clear();
	colorScheme.push_back( GREEN );
	colorScheme.push_back( ( alternative ) ? PURPLE : RED );
	colorScheme.push_back( ( alternative ) ? RED : BLACK );
}

///////////////////////////////////////////////////////////////////////////////
// Check whether the predicted and accepted files can be properly compared.
///////////////////////////////////////////////////////////////////////////////
string StructureComparedImageHandler::validateFiles( string predicted, string accepted, int number ) {

	// For the first structure file, get the necessary validation data.
	// If an error occurred, return it.
	RNA* strand1 = new RNA( predicted.c_str(), FILE_CT );
	ErrorChecker<RNA>* check1 = new ErrorChecker<RNA>( strand1 );
	string strand1Error = "";
	int length1 = -1, linkerIndex1 = -1;
	if( check1->returnError() != "" ) { strand1Error = check1->returnError(); }
	else {
		length1 = strand1->GetSequenceLength();
		numPredictedStructures = strand1->GetStructureNumber();
		for( int i = 1; i <= length1; i++ ) {
			if( strand1->GetNucleotide( i ) == 'I' ) {
				linkerIndex1 = i;
				i += length1;
			}
		}
	}
	delete check1;
	delete strand1;
	if( strand1Error != "" ) { return strand1Error; }

	// For the second structure file, get the necessary validation data.
	// If an error occurred, return it.
	RNA* strand2 = new RNA( accepted.c_str(), FILE_CT );
	ErrorChecker<RNA>* check2 = new ErrorChecker<RNA>( strand2 );
	string strand2Error = "";
	int length2 = -1, linkerIndex2 = -1;
	if( check2->returnError() != "" ) { strand2Error = check2->returnError(); }
	else {
		if( strand2->GetStructureNumber() != 1 ) { strand2Error = "The accepted structures CT file must contain exactly one structure."; }
		else {
			length2 = strand2->GetSequenceLength();
			for( int i = 1; i <= length2; i++ ) {
				if( strand2->GetNucleotide( i ) == 'I' ) {
					linkerIndex2 = i;
					i += length2;
				}
			}
		}
	}
	delete check2;
	delete strand2;
	if( strand2Error != "" ) { return strand2Error; }

	// If the lengths of the sequences don't match, return an error.
	// If a structure number was set to be drawn from the predicted structure, but the number is out of range, return an error.
	if( length1 != length2 ) { return "The predicted and accepted structures are not the same length."; }
	if( ( number != -1 ) && !( number >= 1 && number <= numPredictedStructures ) ) { return "The given structure number does not exist in the predicted structures CT file."; }

	// If the location of the structures' bimolecular linkers are not the same, return an error, as necessary.
	if( linkerIndex1 != linkerIndex2 ) {
		if( linkerIndex1 == -1 ) {
			return "Cannot compare these files; accepted structure is bimolecular and predicted structures are not.";
		} else if( linkerIndex2 == -1 ) {
			return "Cannot compare these files; predicted structures are bimolecular and accepted structure is not.";
		} else {
			return "Cannot compare these files; the bimolecular structures don't have components of the same length.";
		}
	}

	// Return an empty string if the files can be compared properly.
	return "";
}
