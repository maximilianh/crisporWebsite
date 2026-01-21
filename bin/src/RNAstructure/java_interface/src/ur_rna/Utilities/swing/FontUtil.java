package ur_rna.Utilities.swing;

/**
 * Font Utilities
 */
import java.awt.*;
import java.awt.font.FontRenderContext;
import java.awt.font.TextLayout;
import java.text.DecimalFormat;
import java.util.*;
import java.util.Queue;

public abstract class FontUtil {
    private static final FontRenderContext FONT_RENDER_CONTEXT = new FontRenderContext(null, false, true);
    private static final String FONT_TEST_STRING = "The quick brown fox jumps over the lazy dog";
    public static final Font defaultFont = Font.decode(null);

    private static class FontExpressivenessComparator implements Comparator<Font> {
        private static final int[] STYLES = {
                Font.PLAIN, Font.ITALIC, Font.BOLD, Font.BOLD | Font.ITALIC
        };
        public int compare(Font font1, Font font2) {
            if (font1 == font2) {
                return 0;
            }
            Set<String> variantNames1 = new HashSet<>();
            Set<String> variantNames2 = new HashSet<>();
            for (int style : STYLES) {
                variantNames1.add(font1.deriveFont(style).getPSName());
                variantNames2.add(font2.deriveFont(style).getPSName());
            }
            if (variantNames1.size() < variantNames2.size()) {
                return 1;
            } else if (variantNames1.size() > variantNames2.size()) {
                return -1;
            }
            return font1.getName().compareTo(font2.getName());
        }
    }

    private static final FontExpressivenessComparator FONT_EXPRESSIVENESS_COMPARATOR =
            new FontExpressivenessComparator();

    private static boolean isLogicalFontFamily(String family) {
        return (Font.DIALOG.equals(family) ||
                Font.DIALOG_INPUT.equals(family) ||
                Font.SANS_SERIF.equals(family) ||
                Font.SERIF.equals(family) ||
                Font.MONOSPACED.equals(family));
    }

    /**
     * Try to guess physical font from the properties of a logical font, like
     * "Dialog", "Serif", "Monospaced" etc.
     * @param logicalFont Logical font object.
     * @param testText Text used to determine font properties.
     * @return An object of the first matching physical font. The original font
     * object is returned if it was a physical font or no font matched.
     */
    public static Font getPhysicalFont(Font logicalFont, String testText) {
        String logicalFamily = logicalFont.getFamily();
        if (!isLogicalFontFamily(logicalFamily))
            return logicalFont;

        final TextLayout logicalLayout = new TextLayout(testText, logicalFont, FONT_RENDER_CONTEXT);

        // Create a list of matches sorted by font expressiveness (in descending order)
        Queue<Font> physicalFonts = new PriorityQueue<Font>(1, FONT_EXPRESSIVENESS_COMPARATOR);

        Font[] allPhysicalFonts = GraphicsEnvironment.getLocalGraphicsEnvironment().getAllFonts();
        for (Font physicalFont : allPhysicalFonts) {
            String physicalFamily = physicalFont.getFamily();
            // Skip logical fonts
            if (isLogicalFontFamily(physicalFamily)) {
                continue;
            }

            // Derive identical variant of physical font
            physicalFont = physicalFont.deriveFont(
                    logicalFont.getStyle(), logicalFont.getSize2D());
            TextLayout physicalLayout =
                    new TextLayout(testText, physicalFont, FONT_RENDER_CONTEXT);

            // Compare various properties of physical and logical font
            if (physicalLayout.getBounds().equals(logicalLayout.getBounds()) &&
                    physicalLayout.getAscent() == logicalLayout.getAscent() &&
                    physicalLayout.getDescent() == logicalLayout.getDescent() &&
                    physicalLayout.getLeading() == logicalLayout.getLeading() &&
                    physicalLayout.getAdvance() == logicalLayout.getAdvance() &&
                    physicalLayout.getVisibleAdvance() == logicalLayout.getVisibleAdvance()) {
                // Store matching font in list
                physicalFonts.add(physicalFont);
            }
        }

        // Return a valid font even when no matching font could be found
        if (physicalFonts.isEmpty()) {
            return logicalFont;
        }

        return physicalFonts.poll();
    }

    public static Font getPhysicalFont(Font logicalFont) {
        return getPhysicalFont(logicalFont, FONT_TEST_STRING);
    }

    private static DecimalFormat _fontSizeFormatter = new DecimalFormat("#.##");
    public static String encode(Font f) {
        String  strStyle;
        if (f.isBold()) {
            strStyle = f.isItalic() ? "bolditalic" : "bold";
        } else {
            strStyle = f.isItalic() ? "italic" : "plain";
        }

        float size = f.getSize2D();
        return f.getFamily() + '-' + strStyle + '-' + ((size == (int)size) ? Integer.toString((int)size) :  _fontSizeFormatter.format(f));
    }
}
