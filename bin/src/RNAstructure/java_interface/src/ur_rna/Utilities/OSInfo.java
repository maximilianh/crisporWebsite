/*
 * (c) 2012 Mathews Lab, University of Rochester Medical Center.
 *
 * This software is part of a group specifically designed for the RNAstructure
 * secondary structure prediction program and its related applications.
 */

package ur_rna.Utilities;

import javax.swing.*;
import java.io.File;
import java.io.Serializable;

/**
 * A class that handles checking whether an OS is a Mac, and setting general
 * properties if that is the case.
 *
 * @author Richard Watson  (updated detection)
 */
public class OSInfo
	implements Serializable {
	private static final long serialVersionUID = 20120802;
	private static final OSType osType;

	static {
		osType = getOSType(System.getProperty("os.name"));
	}

	public static String getDocumentsDirPath() {
		File f = getDocumentsDir();
		return f == null ? null : f.getPath();
	}
	public static File getDocumentsDir() {
		File home = new File(System.getProperty("user.home"));
		File found = PathTools.getFirstExisting(home,true, false, "Documents", "My Documents", "Docs", "docs");
		return found == null ? home : found;
	}

	public enum OSType {
		WINDOWS,
		LINUX,
		SOLARIS,
		MACOSX,
		UNKNOWN
	}

	private static OSType getOSType(String osName) {
		if(osName != null) {
			osName = osName.toLowerCase();
			if(osName.contains("windows"))
				return OSType.WINDOWS;
			if(osName.contains("linux"))
				return OSType.LINUX;
			if(osName.contains("solaris") || osName.contains("sunos"))
				return OSType.SOLARIS;
			if(osName.contains("os x"))
				return OSType.MACOSX;
		}
		return OSType.UNKNOWN;
	}

	/**
	 * Check the OS to see if it's a Mac, and if it is, set a few system
	 * properties to make things more "Mac-like."
	 */
	public static void applyNativeLookAndFeel() {
		if( isMac() ) {
			System.setProperty("com.apple.mrj.application.live-resize", "true");
		}
		// Set the look and feel for the current OS.
		String lookAndFeel = UIManager.getSystemLookAndFeelClassName();
		try {
			UIManager.setLookAndFeel( lookAndFeel );
		} catch (Exception ex) {
			System.err.println("Error setting native Look-And-Feel to " + lookAndFeel);
			ex.printStackTrace();
		}
	}
	public static void useNativeMenus() {
		if( isMac() ) {
			System.setProperty("apple.laf.useScreenMenuBar", "true");
		}
	}

	public static OSType getOSType() { return osType; }
	/**
	 * Get whether the current OS is a Mac.
	 * @return   True if the OS is a Mac, false if not.
	 */
	public static boolean isMac() {
		//return ( System.getProperty( "mrj.version" ) != null );
		return osType == OSType.MACOSX;
	}

	public static boolean isWin() {
		return osType == OSType.WINDOWS;
	}

	public static boolean isNix() {
		return osType == OSType.LINUX;
	}
}
