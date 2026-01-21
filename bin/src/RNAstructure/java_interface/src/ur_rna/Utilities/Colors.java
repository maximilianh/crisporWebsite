package ur_rna.Utilities;

import java.awt.*;
import java.lang.reflect.Field;
import java.util.Map;
import java.util.TreeMap;

/**
 * List of common Colors and color utilities.
 */
public class Colors {
    public static final Color AliceBlue = new Color(0xF0, 0xF8, 0xFF);
    public static final Color AntiqueWhite = new Color(0xFA, 0xEB, 0xD7);
    public static final Color Aqua = new Color(0x00, 0xFF, 0xFF);
    public static final Color Aquamarine = new Color(0x7F, 0xFF, 0xD4);
    public static final Color Azure = new Color(0xF0, 0xFF, 0xFF);
    public static final Color Beige = new Color(0xF5, 0xF5, 0xDC);
    public static final Color Bisque = new Color(0xFF, 0xE4, 0xC4);
    public static final Color Black = new Color(0x00, 0x00, 0x00);
    public static final Color BlanchedAlmond = new Color(0xFF, 0xEB, 0xCD);
    public static final Color Blue = new Color(0x00, 0x00, 0xFF);
    public static final Color BlueViolet = new Color(0x8A, 0x2B, 0xE2);
    public static final Color Brown = new Color(0xA5, 0x2A, 0x2A);
    public static final Color BurlyWood = new Color(0xDE, 0xB8, 0x87);
    public static final Color CadetBlue = new Color(0x5F, 0x9E, 0xA0);
    public static final Color Chartreuse = new Color(0x7F, 0xFF, 0x00);
    public static final Color Chocolate = new Color(0xD2, 0x69, 0x1E);
    public static final Color Coral = new Color(0xFF, 0x7F, 0x50);
    public static final Color CornflowerBlue = new Color(0x64, 0x95, 0xED);
    public static final Color Cornsilk = new Color(0xFF, 0xF8, 0xDC);
    public static final Color Crimson = new Color(0xDC, 0x14, 0x3C);
    public static final Color Cyan = new Color(0x00, 0xFF, 0xFF);
    public static final Color DarkBlue = new Color(0x00, 0x00, 0x8B);
    public static final Color DarkCyan = new Color(0x00, 0x8B, 0x8B);
    public static final Color DarkGoldenRod = new Color(0xB8, 0x86, 0x0B);
    public static final Color DarkGray = new Color(0xA9, 0xA9, 0xA9);
    public static final Color DarkGreen = new Color(0x00, 0x64, 0x00);
    public static final Color DarkKhaki = new Color(0xBD, 0xB7, 0x6B);
    public static final Color DarkMagenta = new Color(0x8B, 0x00, 0x8B);
    public static final Color DarkOliveGreen = new Color(0x55, 0x6B, 0x2F);
    public static final Color DarkOrange = new Color(0xFF, 0x8C, 0x00);
    public static final Color DarkOrchid = new Color(0x99, 0x32, 0xCC);
    public static final Color DarkRed = new Color(0x8B, 0x00, 0x00);
    public static final Color DarkSalmon = new Color(0xE9, 0x96, 0x7A);
    public static final Color DarkSeaGreen = new Color(0x8F, 0xBC, 0x8F);
    public static final Color DarkSlateBlue = new Color(0x48, 0x3D, 0x8B);
    public static final Color DarkSlateGray = new Color(0x2F, 0x4F, 0x4F);
    public static final Color DarkTurquoise = new Color(0x00, 0xCE, 0xD1);
    public static final Color DarkViolet = new Color(0x94, 0x00, 0xD3);
    public static final Color DeepPink = new Color(0xFF, 0x14, 0x93);
    public static final Color DeepSkyBlue = new Color(0x00, 0xBF, 0xFF);
    public static final Color DimGray = new Color(0x69, 0x69, 0x69);
    public static final Color DodgerBlue = new Color(0x1E, 0x90, 0xFF);
    public static final Color FireBrick = new Color(0xB2, 0x22, 0x22);
    public static final Color FloralWhite = new Color(0xFF, 0xFA, 0xF0);
    public static final Color ForestGreen = new Color(0x22, 0x8B, 0x22);
    public static final Color Fuchsia = new Color(0xFF, 0x00, 0xFF);
    public static final Color Gainsboro = new Color(0xDC, 0xDC, 0xDC);
    public static final Color GhostWhite = new Color(0xF8, 0xF8, 0xFF);
    public static final Color Gold = new Color(0xFF, 0xD7, 0x00);
    public static final Color GoldenRod = new Color(0xDA, 0xA5, 0x20);
    public static final Color Gray = new Color(0x80, 0x80, 0x80);
    public static final Color Green = new Color(0x00, 0x80, 0x00);
    public static final Color GreenYellow = new Color(0xAD, 0xFF, 0x2F);
    public static final Color HoneyDew = new Color(0xF0, 0xFF, 0xF0);
    public static final Color HotPink = new Color(0xFF, 0x69, 0xB4);
    public static final Color IndianRed = new Color(0xCD, 0x5C, 0x5C);
    public static final Color Indigo = new Color(0x4B, 0x00, 0x82);
    public static final Color Ivory = new Color(0xFF, 0xFF, 0xF0);
    public static final Color Khaki = new Color(0xF0, 0xE6, 0x8C);
    public static final Color Lavender = new Color(0xE6, 0xE6, 0xFA);
    public static final Color LavenderBlush = new Color(0xFF, 0xF0, 0xF5);
    public static final Color LawnGreen = new Color(0x7C, 0xFC, 0x00);
    public static final Color LemonChiffon = new Color(0xFF, 0xFA, 0xCD);
    public static final Color LightBlue = new Color(0xAD, 0xD8, 0xE6);
    public static final Color LightCoral = new Color(0xF0, 0x80, 0x80);
    public static final Color LightCyan = new Color(0xE0, 0xFF, 0xFF);
    public static final Color LightGoldenRodYellow = new Color(0xFA, 0xFA, 0xD2);
    public static final Color LightGray = new Color(0xD3, 0xD3, 0xD3);
    public static final Color LightGreen = new Color(0x90, 0xEE, 0x90);
    public static final Color LightPink = new Color(0xFF, 0xB6, 0xC1);
    public static final Color LightSalmon = new Color(0xFF, 0xA0, 0x7A);
    public static final Color LightSeaGreen = new Color(0x20, 0xB2, 0xAA);
    public static final Color LightSkyBlue = new Color(0x87, 0xCE, 0xFA);
    public static final Color LightSlateGray = new Color(0x77, 0x88, 0x99);
    public static final Color LightSteelBlue = new Color(0xB0, 0xC4, 0xDE);
    public static final Color LightYellow = new Color(0xFF, 0xFF, 0xE0);
    public static final Color Lime = new Color(0x00, 0xFF, 0x00);
    public static final Color LimeGreen = new Color(0x32, 0xCD, 0x32);
    public static final Color Linen = new Color(0xFA, 0xF0, 0xE6);
    public static final Color Magenta = new Color(0xFF, 0x00, 0xFF);
    public static final Color Maroon = new Color(0x80, 0x00, 0x00);
    public static final Color MediumAquaMarine = new Color(0x66, 0xCD, 0xAA);
    public static final Color MediumBlue = new Color(0x00, 0x00, 0xCD);
    public static final Color MediumOrchid = new Color(0xBA, 0x55, 0xD3);
    public static final Color MediumPurple = new Color(0x93, 0x70, 0xDB);
    public static final Color MediumSeaGreen = new Color(0x3C, 0xB3, 0x71);
    public static final Color MediumSlateBlue = new Color(0x7B, 0x68, 0xEE);
    public static final Color MediumSpringGreen = new Color(0x00, 0xFA, 0x9A);
    public static final Color MediumTurquoise = new Color(0x48, 0xD1, 0xCC);
    public static final Color MediumVioletRed = new Color(0xC7, 0x15, 0x85);
    public static final Color MidnightBlue = new Color(0x19, 0x19, 0x70);
    public static final Color MintCream = new Color(0xF5, 0xFF, 0xFA);
    public static final Color MistyRose = new Color(0xFF, 0xE4, 0xE1);
    public static final Color Moccasin = new Color(0xFF, 0xE4, 0xB5);
    public static final Color NavajoWhite = new Color(0xFF, 0xDE, 0xAD);
    public static final Color Navy = new Color(0x00, 0x00, 0x80);
    public static final Color OldLace = new Color(0xFD, 0xF5, 0xE6);
    public static final Color Olive = new Color(0x80, 0x80, 0x00);
    public static final Color OliveDrab = new Color(0x6B, 0x8E, 0x23);
    public static final Color Orange = new Color(0xFF, 0xA5, 0x00);
    public static final Color OrangeRed = new Color(0xFF, 0x45, 0x00);
    public static final Color Orchid = new Color(0xDA, 0x70, 0xD6);
    public static final Color PaleGoldenRod = new Color(0xEE, 0xE8, 0xAA);
    public static final Color PaleGreen = new Color(0x98, 0xFB, 0x98);
    public static final Color PaleTurquoise = new Color(0xAF, 0xEE, 0xEE);
    public static final Color PaleVioletRed = new Color(0xDB, 0x70, 0x93);
    public static final Color PapayaWhip = new Color(0xFF, 0xEF, 0xD5);
    public static final Color PeachPuff = new Color(0xFF, 0xDA, 0xB9);
    public static final Color Peru = new Color(0xCD, 0x85, 0x3F);
    public static final Color Pink = new Color(0xFF, 0xC0, 0xCB);
    public static final Color Plum = new Color(0xDD, 0xA0, 0xDD);
    public static final Color PowderBlue = new Color(0xB0, 0xE0, 0xE6);
    public static final Color Purple = new Color(0x80, 0x00, 0x80);
    public static final Color Red = new Color(0xFF, 0x00, 0x00);
    public static final Color RosyBrown = new Color(0xBC, 0x8F, 0x8F);
    public static final Color RoyalBlue = new Color(0x41, 0x69, 0xE1);
    public static final Color SaddleBrown = new Color(0x8B, 0x45, 0x13);
    public static final Color Salmon = new Color(0xFA, 0x80, 0x72);
    public static final Color SandyBrown = new Color(0xF4, 0xA4, 0x60);
    public static final Color SeaGreen = new Color(0x2E, 0x8B, 0x57);
    public static final Color SeaShell = new Color(0xFF, 0xF5, 0xEE);
    public static final Color Sienna = new Color(0xA0, 0x52, 0x2D);
    public static final Color Silver = new Color(0xC0, 0xC0, 0xC0);
    public static final Color SkyBlue = new Color(0x87, 0xCE, 0xEB);
    public static final Color SlateBlue = new Color(0x6A, 0x5A, 0xCD);
    public static final Color SlateGray = new Color(0x70, 0x80, 0x90);
    public static final Color Snow = new Color(0xFF, 0xFA, 0xFA);
    public static final Color SpringGreen = new Color(0x00, 0xFF, 0x7F);
    public static final Color SteelBlue = new Color(0x46, 0x82, 0xB4);
    public static final Color Tan = new Color(0xD2, 0xB4, 0x8C);
    public static final Color Teal = new Color(0x00, 0x80, 0x80);
    public static final Color Thistle = new Color(0xD8, 0xBF, 0xD8);
    public static final Color Tomato = new Color(0xFF, 0x63, 0x47);
    public static final Color Turquoise = new Color(0x40, 0xE0, 0xD0);
    public static final Color Violet = new Color(0xEE, 0x82, 0xEE);
    public static final Color Wheat = new Color(0xF5, 0xDE, 0xB3);
    public static final Color White = new Color(0xFF, 0xFF, 0xFF);
    public static final Color WhiteSmoke = new Color(0xF5, 0xF5, 0xF5);
    public static final Color Yellow = new Color(0xFF, 0xFF, 0x00);
    public static final Color YellowGreen = new Color(0x9A, 0xCD, 0x32);

    private static Map<String,Color> _nameToColor;
    public static String getName(Color c) {
        int rgb = c.getRGB();
        if (_nameToColor == null) buildColors();
        for(Map.Entry<String, Color> kv : _nameToColor.entrySet())
            if (kv.getValue().getRGB() == rgb) return kv.getKey();
        // return the hex name
        return '#'+toHex(c);
    }
    public static String toHex(Color c) {return String.format("%02X%02X%02X", c.getRed(), c.getGreen(), c.getBlue()); }
    public static String toHex(Color c, String prefix) {
        return prefix+toHex(c);
    }

    public static Color getColor(String name) { return getColor(name, false); }
    public static Color getColor(String name, Color defaultIfNotFound) {
        Color c = getColor(name, false);
        return c == null ? defaultIfNotFound : c;
    }
    public static Color getColor(String name, boolean decodeIfNotFound) {
        name = name.trim();
        if (name.isEmpty() || name.charAt(0)=='#')
            return decodeHexRGB(name);
        if (_nameToColor == null) buildColors();
        Color c = _nameToColor.get(name);
        if (c == null && decodeIfNotFound)
            return decodeHexRGB(name);
        return c;
    }

    /**
     * This is the same as {@link Color#decode(String)} except that
     * it does NOT throw a {@link NumberFormatException}, but instead returns null if decoding fails.
     * @param rgb An HTML color string with or without a # or 0x prefix, e.g.: #FF00FF, FF00FF, or 0xFF00FF
     * @return A color if the rgb value is valid or null otherwise.
     */
    public static Color decodeHexRGB(String rgb) {
        try {
            if (rgb == null || rgb.isEmpty()) return null;
            if (rgb.length() == 1 || rgb.charAt(0) != '#' && !(rgb.charAt(0)=='0' && (rgb.charAt(1)=='x'||rgb.charAt(1)=='X')))
                rgb = "#"+rgb; // prepend radix indicator
            return new Color(Integer.decode(rgb));
        } catch (NumberFormatException ex) {
            return null;
        }
    }

    private static void buildColors() {
        _nameToColor = new TreeMap<>(String.CASE_INSENSITIVE_ORDER);
        for (Field f : Colors.class.getFields())
            if (f.getType() == Color.class)
                try {
                    _nameToColor.put(f.getName(), (Color) f.get(null));
                } catch (IllegalAccessException ex) {
                    ex.printStackTrace();
                }
    }

    public static Color setAlpha(Color c, int alpha0to255) { return new Color(alpha0to255 << 24 | (0x00FFFFFF & c.getRGB()), true); }
    public static Color setAlpha(Color c, float alpha0to1) { return setAlpha(c, Math.round(alpha0to1*255)); }
}
