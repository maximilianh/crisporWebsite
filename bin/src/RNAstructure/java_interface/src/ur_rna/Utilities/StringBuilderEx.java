package ur_rna.Utilities;

import java.io.Serializable;
import java.util.stream.IntStream;

/**
 * Delegates most functionality to the {@link StringBuilder} class, while adding a few additional methods.
 */
public class StringBuilderEx implements Serializable, Appendable, CharSequence {
    static final long serialVersionUID = 42L;

    private final StringBuilder _sb;

    /**
     * Use the given StringBuilder as the internal builder.
     * This is useful if a StringBuilder has already been created.
     * @param innerStringBuilder
     */
    public StringBuilderEx(final StringBuilder innerStringBuilder) {
        this._sb = innerStringBuilder;
    }

    /**
     * Get the internal StringBuilder.
     * It is safe to call any methods on this, because {@link StringBuilderEx} does not store any state data
     * that would be altered by normal StringBuilder use.
     *
     * @return The internal StringBuilder.
     */
    public StringBuilder getStringBuilder() {
        return _sb;
    }

    /**
     * Remove the specified number of characters from the end of the string.
     * This method is lenient in that no error is thrown when the StringBuilder's length
     * is less than the number of characters specified.
     *
     * @param numChars The number of characters to remove from the end.
     * @return This {@link StringBuilderEx} (to provide a fluent interface).
     */
    public StringBuilderEx truncate(int numChars) {
        int len = _sb.length() - numChars;
        _sb.setLength(Math.max(len, 0));
        return this;
    }

    /**
     * Remove one or more instances of the specified character from the end of the StringBuilder.
     * @param remove The character to remove, if it is found at the end of the StringBuilder.
     * @return This {@link StringBuilderEx} (to provide a fluent interface).
     */
    public StringBuilderEx trimEnd(char remove) { return trimEnd(remove, -1); }
    /**
     * Remove one or more instances of the specified character from the end of the StringBuilder.
     * @param remove The character to remove, if it is found at the end of the StringBuilder.
     * @param limit The maximum number of characters to remove from the end.
     *              If this is 0, no trimming will be performed.
     *              A value of -1 indicates no limit.
     * @return This {@link StringBuilderEx} (to provide a fluent interface).
     *
     */
    public StringBuilderEx trimEnd(final char remove, int limit) {
        int pos = _sb.length() - 1;
        if (limit == 0) return this;

        // the limit argument specifies the maximum number of characters to remove,
        //   but for convenience this is converted into the 1 less than the MINIMUM returnable length.
        if (limit != -1)
            limit = pos - limit; // i.e. continue the loop until pos == length() - 1 - limit

        while(pos != limit && _sb.charAt(pos) == remove)
            pos--;
        pos++; //back to length
        if (pos != _sb.length())
            _sb.setLength(pos);
        return this;
    }

    /**
     * Remove one or more instances of any of the specified characters from the end of the StringBuilder.
     * @param remove A list of characters to remove, if it any are found at the end of the StringBuilder.
     * @return This {@link StringBuilderEx} (to provide a fluent interface).
     */
    public StringBuilderEx trimEnd(char... remove) {
        return trimEnd(remove, -1);
    }
    /**
     * Remove one or more instances of any of the specified characters from the end of the StringBuilder.
     * @param remove A list of characters to remove, if it any are found at the end of the StringBuilder.
     * @return This {@link StringBuilderEx} (to provide a fluent interface).
     * @param limit The maximum number of characters to remove from the end.
     *              If this is 0, no trimming will be performed.
     *              A value of -1 indicates no limit.
     */
    public StringBuilderEx trimEnd(char[] remove, int limit) {
        int pos = _sb.length();
        int remLen = remove.length;
        if (remLen == 0 || limit==0) return this;

        // the limit argument specifies the maximum number of characters to remove,
        //   but for convenience this is converted into the 1 less than the MINIMUM returnable length.
        if (limit != -1)
            limit = pos - limit; // i.e. continue the loop until pos == length() - 1 - limit

        search:
        while(--pos != limit) {
            char c = _sb.charAt(pos);
            for (int i = remLen; i >= 0; i--)
                if (c == remove[i])
                    continue search;
            // if we didn't continue above, then no character was found, so exit the search
            break;
        }

        pos++; //back to length
        if (pos != _sb.length())
            _sb.setLength(pos);
        return this;
    }

    public StringBuilderEx appendIfEmpty(final String str) { if (_sb.length()==0) _sb.append(str); return this;}
    public StringBuilderEx appendIfNotEmpty(final String str) { if (_sb.length()>0) _sb.append(str); return this;}

    /**
     * Appends the given separator, but only if there is existing text and
     * it does not already end with the separator.
     *
     * @param separator The separator string to append.
     * @return This {@link StringBuilderEx} (to provide a fluent interface).
     */
    public StringBuilderEx appendSeparator(final String separator) {
        if (_sb.length() != 0 && !endsWith(separator))
            _sb.append(separator);
        return this;
    }
//    /**
//     * Appends the given separator, but only if there is existing text and
//     * it does not already end with the separator.
//     *
//     * @param separator The separator string to append.
//     * @return This {@link StringBuilderEx} (to provide a fluent interface).
//     */
//    public StringBuilderEx appendAfterSeparator(final String append, final String separator) {
//        if (_sb.length() != 0 && !endsWith(separator))
//            _sb.append(separator);
//        _sb.append(append);
//        return this;
//    }
    public boolean startsWith(String test) {
        return endsWith(test, false);
    }
    public boolean startsWith(String test, boolean ignoreCase) {
        int len = _sb.length();
        int end = test.length() - 1;
        return len > end &&
                (ignoreCase ?
                        _sb.substring(0, end).equalsIgnoreCase(test) :
                        _sb.substring(0, end).equals(test)
                );
    }
    public boolean endsWith(String test) {
        return endsWith(test, false);
    }
    public boolean endsWith(String test, boolean ignoreCase) {
        int len = _sb.length();
        int start = len - test.length();
        return start >= 0 &&
                (ignoreCase ?
                    _sb.substring(start, len -1).equalsIgnoreCase(test) :
                    _sb.substring(start, len -1).equals(test)
                );
    }

    /* ------------------------------------------------------------------------- */
    /* ---------------- "Delegated" StringBulder Constructors ------------------ */
    /* ------------------------------------------------------------------------- */
    public StringBuilderEx() {
        this._sb = new StringBuilder();
    }

    public StringBuilderEx(final int capacity) {
        this._sb = new StringBuilder(capacity);
    }
    public StringBuilderEx(final String str) {
        this._sb = new StringBuilder(str);
    }
    public StringBuilderEx(final CharSequence seq) {
        this._sb = new StringBuilder(seq);
    }

    /* ------------------------------------------------------------------------- */
    /* ----------------- Delegation of StringBulder Methods -------------------- */
    /* ------------------------------------------------------------------------- */
    public StringBuilderEx append(final float f) {_sb.append(f); return this;}
    public StringBuilderEx append(final Object obj) {_sb.append(obj); return this;}
    public StringBuilderEx append(final String str) {_sb.append(str); return this;}
    /**
     * Appends the specified {@code StringBuffer} to this sequence.
     * <p>
     * The characters of the {@code StringBuffer} argument are appended,
     * in order, to this sequence, increasing the
     * length of this sequence by the length of the argument.
     * If {@code _sb} is {@code null}, then the four characters
     * {@code "null"} are appended to this sequence.
     * <p>
     * Let <i>n</i> be the length of this character sequence just prior to
     * execution of the {@code append} method. Then the character at index
     * <i>k</i> in the new character sequence is equal to the character at
     * index <i>k</i> in the old character sequence, if <i>k</i> is less than
     * <i>n</i>; otherwise, it is equal to the character at index <i>k-n</i>
     * in the argument {@code _sb}.
     *
     * @param   sb   the {@code StringBuffer} to append.
     * @return a reference to this object.
     */
    public StringBuilderEx append(final StringBuffer sb) {_sb.append(sb); return this;}
    @Override
    public StringBuilderEx append(final CharSequence s) {_sb.append(s); return this;}
    /**
     * @throws IndexOutOfBoundsException {@inheritDoc}
     * @param s
     * @param start
     * @param end
     */
    @Override
    public StringBuilderEx append(final CharSequence s, final int start, final int end) {_sb.append(s, start, end); return this;}
    public StringBuilderEx append(final char[] str) {_sb.append(str); return this;}
    /**
     * @throws IndexOutOfBoundsException {@inheritDoc}
     * @param str
     * @param offset
     * @param len
     */
    public StringBuilderEx append(final char[] str, final int offset, final int len) {_sb.append(str, offset, len); return this;}
    public StringBuilderEx append(final boolean b) {_sb.append(b); return this;}
    @Override
    public StringBuilderEx append(final char c) {_sb.append(c); return this;}
    public StringBuilderEx append(final int i) {_sb.append(i); return this;}
    public StringBuilderEx append(final long lng) {_sb.append(lng); return this;}
    public StringBuilderEx append(final double d) {_sb.append(d); return this;}
    /**
     * @since 1.5
     * @param codePoint
     */
    public StringBuilderEx appendCodePoint(final int codePoint) {_sb.appendCodePoint(codePoint); return this;}
    /**
     * @throws StringIndexOutOfBoundsException {@inheritDoc}
     * @param start
     * @param end
     */
    public StringBuilderEx delete(final int start, final int end) {_sb.delete(start, end); return this;}
    /**
     * @throws StringIndexOutOfBoundsException {@inheritDoc}
     * @param index
     */
    public StringBuilderEx deleteCharAt(final int index) {_sb.deleteCharAt(index); return this;}
    /**
     * @throws StringIndexOutOfBoundsException {@inheritDoc}
     * @param start
     * @param end
     * @param str
     */
    public StringBuilderEx replace(final int start, final int end, final String str) {_sb.replace(start, end, str); return this;}
    /**
     * @throws StringIndexOutOfBoundsException {@inheritDoc}
     * @param index
     * @param str
     * @param offset
     * @param len
     */
    public StringBuilderEx insert(final int index, final char[] str, final int offset, final int len) {_sb.insert(index, str, offset, len); return this;}
    /**
     * @throws StringIndexOutOfBoundsException {@inheritDoc}
     * @param offset
     * @param obj
     */
    public StringBuilderEx insert(final int offset, final Object obj) {_sb.insert(offset, obj); return this;}
    /**
     * @throws StringIndexOutOfBoundsException {@inheritDoc}
     * @param offset
     * @param str
     */
    public StringBuilderEx insert(final int offset, final String str) {_sb.insert(offset, str); return this;}
    /**
     * @throws StringIndexOutOfBoundsException {@inheritDoc}
     * @param offset
     * @param str
     */
    public StringBuilderEx insert(final int offset, final char[] str) {_sb.insert(offset, str); return this;}
    /**
     * @throws IndexOutOfBoundsException {@inheritDoc}
     * @param dstOffset
     * @param s
     */
    public StringBuilderEx insert(final int dstOffset, final CharSequence s) {_sb.insert(dstOffset, s); return this;}
    /**
     * @throws IndexOutOfBoundsException {@inheritDoc}
     * @param dstOffset
     * @param s
     * @param start
     * @param end
     */
    public StringBuilderEx insert(final int dstOffset, final CharSequence s, final int start, final int end) {_sb.insert(dstOffset, s, start, end); return this;}
    /**
     * @throws StringIndexOutOfBoundsException {@inheritDoc}
     * @param offset
     * @param b
     */
    public StringBuilderEx insert(final int offset, final boolean b) {_sb.insert(offset, b); return this;}
    /**
     * @throws IndexOutOfBoundsException {@inheritDoc}
     * @param offset
     * @param c
     */
    public StringBuilderEx insert(final int offset, final char c) {_sb.insert(offset, c); return this;}
    /**
     * @throws StringIndexOutOfBoundsException {@inheritDoc}
     * @param offset
     * @param i
     */
    public StringBuilderEx insert(final int offset, final int i) {_sb.insert(offset, i); return this;}
    /**
     * @throws StringIndexOutOfBoundsException {@inheritDoc}
     * @param offset
     * @param l
     */
    public StringBuilderEx insert(final int offset, final long l) {_sb.insert(offset, l); return this;}
    /**
     * @throws StringIndexOutOfBoundsException {@inheritDoc}
     * @param offset
     * @param f
     */
    public StringBuilderEx insert(final int offset, final float f) {_sb.insert(offset, f); return this;}
    /**
     * @throws StringIndexOutOfBoundsException {@inheritDoc}
     * @param offset
     * @param d
     */
    public StringBuilderEx insert(final int offset, final double d) {_sb.insert(offset, d); return this;}
    public int indexOf(final String str) {return _sb.indexOf(str);}
    public int indexOf(final String str, final int fromIndex) {return _sb.indexOf(str, fromIndex);}
    public int lastIndexOf(final String str) {return _sb.lastIndexOf(str);}
    public int lastIndexOf(final String str, final int fromIndex) {return _sb.lastIndexOf(str, fromIndex);}
    public StringBuilderEx reverse() {_sb.reverse(); return this;}
    @Override
    public String toString() {return _sb.toString();}
    /**
     * Returns the length (character count).
     *
     * @return the length of the sequence of characters currently
     *          represented by this object
     */
    @Override
    public int length() {return _sb.length();}
    /**
     * Returns the current capacity. The capacity is the amount of storage
     * available for newly inserted characters, beyond which an allocation
     * will occur.
     *
     * @return the current capacity
     */
    public int capacity() {return _sb.capacity();}
    /**
     * Ensures that the capacity is at least equal to the specified minimum.
     * If the current capacity is less than the argument, then a new internal
     * array is allocated with greater capacity. The new capacity is the
     * larger of:
     * <ul>
     * <li>The {@code minimumCapacity} argument.
     * <li>Twice the old capacity, plus {@code 2}.
     * </ul>
     * If the {@code minimumCapacity} argument is nonpositive, this
     * method takes no action and simply returns.
     * Note that subsequent operations on this object can reduce the
     * actual capacity below that requested here.
     *
     * @param   minimumCapacity   the minimum desired capacity.
     */
    public void ensureCapacity(final int minimumCapacity) {_sb.ensureCapacity(minimumCapacity);}
    /**
     * Attempts to reduce storage used for the character sequence.
     * If the buffer is larger than necessary to hold its current sequence of
     * characters, then it may be resized to become more space efficient.
     * Calling this method may, but is not required to, affect the value
     * returned by a subsequent call to the {@link #capacity()} method.
     */
    public void trimToSize() {_sb.trimToSize();}
    /**
     * Sets the length of the character sequence.
     * The sequence is changed to a new character sequence
     * whose length is specified by the argument. For every nonnegative
     * index <i>k</i> less than {@code newLength}, the character at
     * index <i>k</i> in the new character sequence is the same as the
     * character at index <i>k</i> in the old sequence if <i>k</i> is less
     * than the length of the old character sequence; otherwise, it is the
     * null character {@code '\u005Cu0000'}.
     *
     * In other words, if the {@code newLength} argument is less than
     * the current length, the length is changed to the specified length.
     * <p>
     * If the {@code newLength} argument is greater than or equal
     * to the current length, sufficient null characters
     * ({@code '\u005Cu0000'}) are appended so that
     * length becomes the {@code newLength} argument.
     * <p>
     * The {@code newLength} argument must be greater than or equal
     * to {@code 0}.
     *
     * @param      newLength   the new length
     * @throws IndexOutOfBoundsException  if the
     *               {@code newLength} argument is negative.
     */
    public void setLength(final int newLength) {_sb.setLength(newLength);}
    /**
     * Returns the {@code char} value in this sequence at the specified index.
     * The first {@code char} value is at index {@code 0}, the next at index
     * {@code 1}, and so on, as in array indexing.
     * <p>
     * The index argument must be greater than or equal to
     * {@code 0}, and less than the length of this sequence.
     *
     * <p>If the {@code char} value specified by the index is a
     * <a href="Character.html#unicode">surrogate</a>, the surrogate
     * value is returned.
     *
     * @param      index   the index of the desired {@code char} value.
     * @return the {@code char} value at the specified index.
     * @throws IndexOutOfBoundsException  if {@code index} is
     *             negative or greater than or equal to {@code length()}.
     */
    @Override
    public char charAt(final int index) {return _sb.charAt(index);}
    /**
     * Returns the character (Unicode code point) at the specified
     * index. The index refers to {@code char} values
     * (Unicode code units) and ranges from {@code 0} to
     * {@link #length()}{@code  - 1}.
     *
     * <p> If the {@code char} value specified at the given index
     * is in the high-surrogate range, the following index is less
     * than the length of this sequence, and the
     * {@code char} value at the following index is in the
     * low-surrogate range, then the supplementary code point
     * corresponding to this surrogate pair is returned. Otherwise,
     * the {@code char} value at the given index is returned.
     *
     * @param      index the index to the {@code char} values
     * @return the code point value of the character at the
     *             {@code index}
     * @exception IndexOutOfBoundsException  if the {@code index}
     *             argument is negative or not less than the length of this
     *             sequence.
     */
    public int codePointAt(final int index) {return _sb.codePointAt(index);}
    /**
     * Returns the character (Unicode code point) before the specified
     * index. The index refers to {@code char} values
     * (Unicode code units) and ranges from {@code 1} to {@link
     * #length()}.
     *
     * <p> If the {@code char} value at {@code (index - 1)}
     * is in the low-surrogate range, {@code (index - 2)} is not
     * negative, and the {@code char} value at {@code (index -
     * 2)} is in the high-surrogate range, then the
     * supplementary code point value of the surrogate pair is
     * returned. If the {@code char} value at {@code index -
     * 1} is an unpaired low-surrogate or a high-surrogate, the
     * surrogate value is returned.
     *
     * @param     index the index following the code point that should be returned
     * @return the Unicode code point value before the given index.
     * @exception IndexOutOfBoundsException if the {@code index}
     *            argument is less than 1 or greater than the length
     *            of this sequence.
     */
    public int codePointBefore(final int index) {return _sb.codePointBefore(index);}
    /**
     * Returns the number of Unicode code points in the specified text
     * range of this sequence. The text range begins at the specified
     * {@code beginIndex} and extends to the {@code char} at
     * index {@code endIndex - 1}. Thus the length (in
     * {@code char}s) of the text range is
     * {@code endIndex-beginIndex}. Unpaired surrogates within
     * this sequence count as one code point each.
     *
     * @param beginIndex the index to the first {@code char} of
     * the text range.
     * @param endIndex the index after the last {@code char} of
     * the text range.
     * @return the number of Unicode code points in the specified text
     * range
     * @exception IndexOutOfBoundsException if the
     * {@code beginIndex} is negative, or {@code endIndex}
     * is larger than the length of this sequence, or
     * {@code beginIndex} is larger than {@code endIndex}.
     */
    public int codePointCount(final int beginIndex, final int endIndex) {return _sb.codePointCount(beginIndex, endIndex);}
    /**
     * Returns the index within this sequence that is offset from the
     * given {@code index} by {@code codePointOffset} code
     * points. Unpaired surrogates within the text range given by
     * {@code index} and {@code codePointOffset} count as
     * one code point each.
     *
     * @param index the index to be offset
     * @param codePointOffset the offset in code points
     * @return the index within this sequence
     * @exception IndexOutOfBoundsException if {@code index}
     *   is negative or larger then the length of this sequence,
     *   or if {@code codePointOffset} is positive and the subsequence
     *   starting with {@code index} has fewer than
     *   {@code codePointOffset} code points,
     *   or if {@code codePointOffset} is negative and the subsequence
     *   before {@code index} has fewer than the absolute value of
     *   {@code codePointOffset} code points.
     */
    public int offsetByCodePoints(final int index, final int codePointOffset) {return _sb.offsetByCodePoints(index, codePointOffset);}
    /**
     * Characters are copied from this sequence into the
     * destination character array {@code dst}. The first character to
     * be copied is at index {@code srcBegin}; the last character to
     * be copied is at index {@code srcEnd-1}. The total number of
     * characters to be copied is {@code srcEnd-srcBegin}. The
     * characters are copied into the subarray of {@code dst} starting
     * at index {@code dstBegin} and ending at index:
     * <pre>{@code
     * dstbegin + (srcEnd-srcBegin) - 1
     * }</pre>
     *
     * @param      srcBegin   start copying at this offset.
     * @param      srcEnd     stop copying at this offset.
     * @param      dst        the array to copy the data into.
     * @param      dstBegin   offset into {@code dst}.
     * @throws IndexOutOfBoundsException  if any of the following is true:
     *             <ul>
     *             <li>{@code srcBegin} is negative
     *             <li>{@code dstBegin} is negative
     *             <li>the {@code srcBegin} argument is greater than
     *             the {@code srcEnd} argument.
     *             <li>{@code srcEnd} is greater than
     *             {@code this.length()}.
     *             <li>{@code dstBegin+srcEnd-srcBegin} is greater than
     *             {@code dst.length}
     *             </ul>
     */
    public void getChars(final int srcBegin, final int srcEnd, final char[] dst, final int dstBegin) {_sb.getChars(srcBegin, srcEnd, dst, dstBegin);}
    /**
     * The character at the specified index is set to {@code ch}. This
     * sequence is altered to represent a new character sequence that is
     * identical to the old character sequence, except that it contains the
     * character {@code ch} at position {@code index}.
     * <p>
     * The index argument must be greater than or equal to
     * {@code 0}, and less than the length of this sequence.
     *
     * @param      index   the index of the character to modify.
     * @param      ch      the new character.
     * @throws IndexOutOfBoundsException  if {@code index} is
     *             negative or greater than or equal to {@code length()}.
     */
    public void setCharAt(final int index, final char ch) {_sb.setCharAt(index, ch);}
    /**
     * Returns a new {@code String} that contains a subsequence of
     * characters currently contained in this character sequence. The
     * substring begins at the specified index and extends to the end of
     * this sequence.
     *
     * @param      start    The beginning index, inclusive.
     * @return The new string.
     * @throws StringIndexOutOfBoundsException  if {@code start} is
     *             less than zero, or greater than the length of this object.
     */
    public String substring(final int start) {return _sb.substring(start);}
    /**
     * Returns a new character sequence that is a subsequence of this sequence.
     *
     * <p> An invocation of this method of the form
     *
     * <pre>{@code
     * _sb.subSequence(begin,&nbsp;end)}</pre>
     *
     * behaves in exactly the same way as the invocation
     *
     * <pre>{@code
     * _sb.substring(begin,&nbsp;end)}</pre>
     *
     * This method is provided so that this class can
     * implement the {@link CharSequence} interface.
     *
     * @param      start   the start index, inclusive.
     * @param      end     the end index, exclusive.
     * @return the specified subsequence.
     *
     * @throws IndexOutOfBoundsException
     *          if {@code start} or {@code end} are negative,
     *          if {@code end} is greater than {@code length()},
     *          or if {@code start} is greater than {@code end}
     */
    @Override
    public CharSequence subSequence(final int start, final int end) {return _sb.subSequence(start, end);}
    /**
     * Returns a new {@code String} that contains a subsequence of
     * characters currently contained in this sequence. The
     * substring begins at the specified {@code start} and
     * extends to the character at index {@code end - 1}.
     *
     * @param      start    The beginning index, inclusive.
     * @param      end      The ending index, exclusive.
     * @return The new string.
     * @throws StringIndexOutOfBoundsException  if {@code start}
     *             or {@code end} are negative or greater than
     *             {@code length()}, or {@code start} is
     *             greater than {@code end}.
     */
    public String substring(final int start, final int end) {return _sb.substring(start, end);}
    /**
     * Returns a stream of {@code int} zero-extending the {@code char} values
     * from this sequence.  Any char which maps to a <a
     * href="{@docRoot}/java/lang/Character.html#unicode">surrogate code
     * point</a> is passed through uninterpreted.
     *
     * <p>If the sequence is mutated while the stream is being read, the
     * result is undefined.
     *
     * @return an IntStream of char values from this sequence
     * @since 1.8
     */
    @Override
    public IntStream chars() {return _sb.chars();}
    /**
     * Returns a stream of code point values from this sequence.  Any surrogate
     * pairs encountered in the sequence are combined as if by {@linkplain
     * Character#toCodePoint Character.toCodePoint} and the result is passed
     * to the stream. Any other code units, including ordinary BMP characters,
     * unpaired surrogates, and undefined code units, are zero-extended to
     * {@code int} values which are then passed to the stream.
     *
     * <p>If the sequence is mutated while the stream is being read, the result
     * is undefined.
     *
     * @return an IntStream of Unicode code points from this sequence
     * @since 1.8
     */
    @Override
    public IntStream codePoints() {return _sb.codePoints();}
}
