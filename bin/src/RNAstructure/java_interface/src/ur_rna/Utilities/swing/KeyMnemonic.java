package ur_rna.Utilities.swing;

import java.awt.event.KeyEvent;

/**
 * @author Richard M. Watson
 */
public final class KeyMnemonic {
//    private final static Pattern mnemonicRegex = Pattern.compile("(?<!&)&(\\w)");
//    private final static Pattern mnemonicEscapeRegex = Pattern.compile("&&");
    public static String stripMnemonics(String title) {
        if (title == null) return null;
        StringBuilder sb = new StringBuilder(title.length());
        // return the string with all double-ampersands replaced with single ampersands.
        // If there are any single-ampersands, remove only the first one.
        boolean prevAmp = false, foundMnemonic = false;
        for (char c : title.toCharArray()) {
            if (c == '&') {
                prevAmp = !prevAmp;
                if (prevAmp)
                    continue; // i.e. if prevAmp was false (i.e. the previous character was not '&'), do NOT append the character.
                    // sb.append('&'); // append the second of two adjacent ampersands.
            } else if (prevAmp) {
                // The previous character was an ampersand, but this character is not.
                if (foundMnemonic) // if we already found the mnemonic, append this ampersand (which was skipped in the previous iteration).
                    sb.append('&');
                else
                    foundMnemonic = true;
                prevAmp = false;
            }
            sb.append(c);
        }
        if (foundMnemonic && prevAmp)
            sb.append('&');
        return sb.toString();
    }

//    private static Map<String, Integer> vkMap;
//    private static Map<String, Integer> buildVKMap() {
//        Map<String, Integer> map = new HashMap<String, Integer>(256);
//        Field[] fields = KeyEvent.class.getDeclaredFields();
//        for (Field field : fields) {
//            String name = field.getName();
//            if (name.startsWith("VK_")) try {
//                map.put(name.substring(3).toUpperCase(), // remove "VK_" prefix
//                        field.getInt(null));
//            } catch (IllegalAccessException ex) {
//                //
//            }
//        }
//        return map;
//    }
    public static int getMnemonic(final String title) {
        // return the first character after non-doubled ampersand.
        // ex: &File => "F"
        // ex: PB && &Jelly => "J"
        if (title == null) return 0;
        boolean prevAmp = false;
        for (char c : title.toCharArray()) {
            if (c == '&')
                prevAmp = !prevAmp;
            else if (prevAmp)
                return getVKCode(c);
        }
        return 0;
    }

    public static int getVKCode(char c) {
        c = Character.toUpperCase(c);
        if (c >='A' && c <= 'Z')
            return KeyEvent.VK_A + c - 'A';
        if (c >='0' && c <= '9')
            return KeyEvent.VK_0 + c - '0';
        return 0;
    }
    /**
     * Returns the index in the string at which the mnemonic character is found.
     * The index represents the position as it would be in the string after processing by {@link #stripMnemonics}
     * @param title
     * @return Returns -1 if no mnemonic was found. Otherwise returns the position of the mnemonic character (as it would be in the text after {@link #stripMnemonics} is applied.)
     */
    public static int getMnemonicIndex(final String title) {
        // return the first character after non-doubled ampersand.
        // ex: &File => "F"
        // ex: PB && &Jelly => "J"
        if (title == null) return -1;
        boolean prevAmp = false;
        int pos = -1;
        for (char c : title.toCharArray()) {
            pos++;
            if (c == '&') {
                if (!prevAmp)
                    pos--; // the second amp it would be removed stripMnemonics
                prevAmp = !prevAmp;
            } else if (prevAmp)
                return pos;
        }
        return -1;
    }
}
