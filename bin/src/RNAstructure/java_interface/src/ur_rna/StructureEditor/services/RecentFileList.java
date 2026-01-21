package ur_rna.StructureEditor.services;

import ur_rna.StructureEditor.FileType;

import java.util.*;
import java.util.function.Consumer;
import java.util.prefs.BackingStoreException;
import java.util.prefs.Preferences;

/**
 * List of recent files that prevents duplicate paths.
 */
public class RecentFileList implements Iterable<RecentFileList.RecentFile> {
    private int maxItemCount = 16;
    private LinkedList<RecentFileList.RecentFile> list = new LinkedList<>();
    public List<RecentFile> getFiles() {
        return Collections.unmodifiableList(list);
    }
    public int getMaxItemCount() {
        return maxItemCount;
    }
    public void setMaxItemCount(final int maxItemCount) {
        this.maxItemCount = maxItemCount;
    }
    public int size() { return list.size(); }
    public void clear() {
        list.clear();
    }
    public void remove(final RecentFile recentFile) {
        list.remove(recentFile);
    }
    public static class RecentFile {
        public RecentFile(final String path) { this(path, null); }
        public RecentFile(final String path, final FileType type) {
            this.path = path;
            this.type = type;
        }
        public final String path;
        public final FileType type;
    }
    public boolean contains(String path) {
        return indexOf(path)!=-1;
    }
    public int indexOf(String path) {
        int i = 0;
        for (RecentFile f : list) {
            if (path.equalsIgnoreCase(f.path))
                return i;
            i++;
        }
        return -1;
    }
    public void remove(String path) {
        int i = indexOf(path);
        if (i != -1) list.remove(i);
    }
    public void add(String path) { add(path, null); }
    public void add(String path, FileType type) {
        remove(path);
        prepareAdd();
        list.addFirst(new RecentFile(path, type));
    }
    // make room to add a new item. (i.e. size should be at least one less than max)
    private void prepareAdd() {
        while (list.size() >= maxItemCount)
            list.removeLast();
    }
    public void addToBottom(String path, FileType type) {
        prepareAdd();
        list.addLast(new RecentFile(path, type));
    }

    public void load(Preferences p) {
        list.clear();
        Set<String> set = new LinkedHashSet<>();
        try {
            String[] keys = p.keys();
            Arrays.sort(keys);
            for (String key : keys)
                set.add(p.get(key, null));
        } catch (BackingStoreException ex){
            ex.printStackTrace();
        }
        for(String s : set) {
            int pos = s.indexOf('\t');
            if (pos == -1)
                list.addLast(new RecentFile(s, parseFileType(null, s)));
            else {
                String path = s.substring(0, pos);
                list.addLast(new RecentFile(path, parseFileType(s.substring(pos+1), path)));
            }
        }
    }
    public static FileType parseFileType(final String type, final String path) {
        return FileType.fromAny(type, null, null, path, null);
    }

    public void save(Preferences p) {
        int i = 0; // first item will be 1
        try {
            p.clear();
        } catch (BackingStoreException ex) {
            ex.printStackTrace();
        }
        for (RecentFile f : list) {
            if (i++ > maxItemCount) return;
            String key = "00" + Integer.toString(i);
            if (key.length() > 3)
                key = key.substring(key.length() - 3);
            p.put(key, f.path + "\t" + (f.type==null?"":f.type.name()));
        }
        try {
            p.flush();
        } catch (BackingStoreException ex) {
            ex.printStackTrace();
        }
    }
    @Override
    public Iterator<RecentFile> iterator() {
        return list.iterator();
    }
    @Override
    public void forEach(final Consumer<? super RecentFile> action) {
        list.forEach(action);
    }
    @Override
    public Spliterator<RecentFile> spliterator() {
        return list.spliterator();
    }
}
