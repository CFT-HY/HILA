/**
 * @file hdf5_file.cpp
 * @brief Encoding of the HDF5 file format subset written by hila::hdf5_file
 *
 * Structures used (HDF5 File Format Specification Version 3.0):
 * - superblock version 2 (48 bytes at offset 0, written last)
 * - version 2 object headers ("OHDR") with a single chunk
 * - "new style" compact groups: Link Info, Group Info and Link messages in the group object
 *   header, so that no B-trees or heaps are needed
 * - datasets: Dataspace (v2), Datatype (v1, fixed-point / IEEE float), Fill Value (v3) and
 *   Data Layout (v3, contiguous) messages
 * - Attribute messages (v3)
 * All numbers are little-endian, offsets and lengths 8 bytes.  Metadata checksums use the
 * Jenkins lookup3 hash, as in the HDF5 library.
 */

#include "plumbing/hdf5_file.h"

#include <cstring>

namespace {

constexpr uint64_t undefined_address = 0xffffffffffffffffULL;
constexpr uint64_t superblock_size = 48;
constexpr uint64_t data_alignment = 4096; // start of each dataset's raw data

// HDF5 header message types
constexpr uint8_t msg_dataspace = 0x01;
constexpr uint8_t msg_link_info = 0x02;
constexpr uint8_t msg_datatype = 0x03;
constexpr uint8_t msg_fill_value = 0x05;
constexpr uint8_t msg_link = 0x06;
constexpr uint8_t msg_layout = 0x08;
constexpr uint8_t msg_group_info = 0x0A;
constexpr uint8_t msg_attribute = 0x0C;

using bytes_t = std::vector<uint8_t>;

void put8(bytes_t &b, uint64_t v) {
    b.push_back((uint8_t)v);
}
void put16(bytes_t &b, uint64_t v) {
    for (int i = 0; i < 2; i++)
        b.push_back((uint8_t)(v >> (8 * i)));
}
void put32(bytes_t &b, uint64_t v) {
    for (int i = 0; i < 4; i++)
        b.push_back((uint8_t)(v >> (8 * i)));
}
void put64(bytes_t &b, uint64_t v) {
    for (int i = 0; i < 8; i++)
        b.push_back((uint8_t)(v >> (8 * i)));
}
void append(bytes_t &b, const bytes_t &a) {
    b.insert(b.end(), a.begin(), a.end());
}

/// Jenkins lookup3 hash ("hashlittle"), as H5_checksum_lookup3 in the HDF5 library
uint32_t lookup3(const uint8_t *k, size_t length, uint32_t initval = 0) {

    auto rot = [](uint32_t x, int s) -> uint32_t { return (x << s) | (x >> (32 - s)); };

    uint32_t a, b, c;
    a = b = c = 0xdeadbeef + (uint32_t)length + initval;

    while (length > 12) {
        a += k[0] + ((uint32_t)k[1] << 8) + ((uint32_t)k[2] << 16) + ((uint32_t)k[3] << 24);
        b += k[4] + ((uint32_t)k[5] << 8) + ((uint32_t)k[6] << 16) + ((uint32_t)k[7] << 24);
        c += k[8] + ((uint32_t)k[9] << 8) + ((uint32_t)k[10] << 16) + ((uint32_t)k[11] << 24);
        // mix
        a -= c;  a ^= rot(c, 4);  c += b;
        b -= a;  b ^= rot(a, 6);  a += c;
        c -= b;  c ^= rot(b, 8);  b += a;
        a -= c;  a ^= rot(c, 16); c += b;
        b -= a;  b ^= rot(a, 19); a += c;
        c -= b;  c ^= rot(b, 4);  b += a;
        length -= 12;
        k += 12;
    }

    switch (length) {
    case 12: c += (uint32_t)k[11] << 24; [[fallthrough]];
    case 11: c += (uint32_t)k[10] << 16; [[fallthrough]];
    case 10: c += (uint32_t)k[9] << 8;   [[fallthrough]];
    case 9:  c += k[8];                  [[fallthrough]];
    case 8:  b += (uint32_t)k[7] << 24;  [[fallthrough]];
    case 7:  b += (uint32_t)k[6] << 16;  [[fallthrough]];
    case 6:  b += (uint32_t)k[5] << 8;   [[fallthrough]];
    case 5:  b += k[4];                  [[fallthrough]];
    case 4:  a += (uint32_t)k[3] << 24;  [[fallthrough]];
    case 3:  a += (uint32_t)k[2] << 16;  [[fallthrough]];
    case 2:  a += (uint32_t)k[1] << 8;   [[fallthrough]];
    case 1:  a += k[0]; break;
    case 0:  return c;
    }
    // final
    c ^= b;  c -= rot(b, 14);
    a ^= c;  a -= rot(c, 11);
    b ^= a;  b -= rot(a, 25);
    c ^= b;  c -= rot(b, 16);
    a ^= c;  a -= rot(c, 4);
    b ^= a;  b -= rot(a, 14);
    c ^= b;  c -= rot(b, 24);
    return c;
}

/// Datatype message for numbers
bytes_t encode_datatype(const hila::hdf5_file::number_type &t) {
    bytes_t b;
    if (t.is_float) {
        const bool dbl = (t.size == 8);
        put8(b, 0x11);             // class 1 (floating point), version 1
        put8(b, 0x20);             // little-endian, no padding, mantissa msb implied
        put8(b, dbl ? 63 : 31);    // sign bit location
        put8(b, 0);
        put32(b, t.size);          // size in bytes
        put16(b, 0);               // bit offset
        put16(b, 8 * t.size);      // bit precision
        put8(b, dbl ? 52 : 23);    // exponent location
        put8(b, dbl ? 11 : 8);     // exponent size
        put8(b, 0);                // mantissa location
        put8(b, dbl ? 52 : 23);    // mantissa size
        put32(b, dbl ? 1023 : 127); // exponent bias
    } else {
        put8(b, 0x10);                  // class 0 (fixed point), version 1
        put8(b, t.is_signed ? 0x08 : 0); // little-endian, signed bit
        put8(b, 0);
        put8(b, 0);
        put32(b, t.size);
        put16(b, 0);          // bit offset
        put16(b, 8 * t.size); // bit precision
    }
    return b;
}

/// Datatype message for a fixed-length, null-padded ASCII string of `size` bytes (as written
/// by h5py for numpy byte strings).  The string is stored without a terminating null: VTK's
/// VTKHDF reader compares all `size` bytes, e.g. of the attribute Type = "ImageData".
bytes_t encode_string_datatype(uint64_t size) {
    bytes_t b;
    put8(b, 0x13); // class 3 (string), version 1
    put8(b, 0x01); // null padded, ASCII
    put8(b, 0);
    put8(b, 0);
    put32(b, size);
    return b;
}

/// Dataspace message, version 2.  scalar: no dims
bytes_t encode_dataspace(const std::vector<uint64_t> &dims, bool scalar) {
    bytes_t b;
    put8(b, 2);                       // version
    put8(b, scalar ? 0 : dims.size()); // rank
    put8(b, 0);                       // flags: no max dims
    put8(b, scalar ? 0 : 1);          // type: scalar / simple
    if (!scalar)
        for (auto d : dims)
            put64(b, d);
    return b;
}

bytes_t encode_attribute(const std::string &name, const bytes_t &dtype, const bytes_t &dspace,
                         const bytes_t &data) {
    bytes_t b;
    put8(b, 3); // version
    put8(b, 0); // flags
    put16(b, name.size() + 1);
    put16(b, dtype.size());
    put16(b, dspace.size());
    put8(b, 0); // name encoding ASCII
    b.insert(b.end(), name.begin(), name.end());
    put8(b, 0);
    append(b, dtype);
    append(b, dspace);
    append(b, data);
    return b;
}

/// Version 2 object header with a single chunk holding the messages
bytes_t object_header(const std::vector<std::pair<uint8_t, bytes_t>> &messages,
                      const std::string &objname) {
    uint64_t chunk_size = 0;
    for (auto &m : messages) {
        if (m.second.size() > 0xffff)
            hila::error("HDF5 output: header message too large in object " + objname +
                        " (attribute too large?)");
        chunk_size += 4 + m.second.size();
    }

    bytes_t b;
    for (char ch : std::string("OHDR"))
        put8(b, ch);
    put8(b, 2);    // version
    put8(b, 0x02); // flags: chunk #0 size has 4 bytes, no times, no attribute phase change
    put32(b, chunk_size);
    for (auto &m : messages) {
        put8(b, m.first);
        put16(b, m.second.size());
        put8(b, 0); // message flags
        append(b, m.second);
    }
    put32(b, lookup3(b.data(), b.size()));
    return b;
}

std::vector<std::string> split_path(const std::string &path) {
    std::vector<std::string> parts;
    size_t start = 0;
    while (start <= path.size()) {
        size_t end = path.find('/', start);
        if (end == std::string::npos)
            end = path.size();
        if (end > start)
            parts.push_back(path.substr(start, end - start));
        start = end + 1;
    }
    return parts;
}

} // namespace

namespace hila {

hdf5_file::hdf5_file(const std::string &fname, hila::io::output_flags flags) : filename(fname) {

    uint16_t one = 1;
    if (*reinterpret_cast<uint8_t *>(&one) != 1)
        hila::error("hdf5_file: only little-endian systems supported");

    fh = hila::io::open_mpi_file(filename, true);
    is_open = true;
    eof = superblock_size;

    object_t root;
    root.is_group = true;
    objects.push_back(root);

    vtk_mode = hila::io::is_vtk_output(flags, filename);
}

hdf5_file::~hdf5_file() {
    if (is_open && !hila::about_to_finish)
        close();
}

uint64_t hdf5_file::reserve(uint64_t bytes) {
    uint64_t addr = (eof + data_alignment - 1) / data_alignment * data_alignment;
    eof = addr + bytes;
    return addr;
}

int hdf5_file::find_or_create_group(const std::string &path) {
    int cur = 0;
    for (const auto &part : split_path(path)) {
        int next = -1;
        for (int c : objects[cur].children) {
            if (objects[c].name == part) {
                next = c;
                break;
            }
        }
        if (next >= 0) {
            if (!objects[next].is_group)
                hila::error("hdf5_file " + filename + ": " + path +
                            " - path component is a dataset, not a group");
        } else {
            object_t g;
            g.is_group = true;
            g.name = part;
            objects.push_back(g);
            next = objects.size() - 1;
            objects[cur].children.push_back(next);
        }
        cur = next;
    }
    return cur;
}

int hdf5_file::add_dataset(const std::string &path, number_type nt,
                           const std::vector<uint64_t> &dims, uint64_t addr, uint64_t bytes) {
    if (!is_open)
        hila::error("hdf5_file " + filename + ": writing dataset " + path + " to a closed file");

    auto parts = split_path(path);
    if (parts.empty())
        hila::error("hdf5_file " + filename + ": empty dataset name");

    std::string name = parts.back();
    int parent = find_or_create_group(path.substr(0, path.rfind('/') == std::string::npos
                                                         ? 0
                                                         : path.rfind('/')));
    for (int c : objects[parent].children)
        if (objects[c].name == name)
            hila::error("hdf5_file " + filename + ": object " + path + " already exists");

    object_t d;
    d.is_group = false;
    d.name = name;
    d.dtype = encode_datatype(nt);
    d.dspace = encode_dataspace(dims, false);
    d.addr = addr;
    d.size = bytes;
    objects.push_back(d);
    int idx = objects.size() - 1;
    objects[parent].children.push_back(idx);
    return idx;
}

/// Add a dataset on a grid of size sizes_x_first (x first) and number of components ncomp.
/// In VTKHDF mode the dataset goes to /VTKHDF/PointData, and the grid (padded to 3
/// dimensions) must be the same for all datasets.
int hdf5_file::add_grid_dataset(const std::string &path, number_type nt,
                                const std::vector<uint64_t> &sizes_x_first,
                                const std::vector<int64_t> &origin_x_first, size_t ncomp,
                                uint64_t addr, uint64_t bytes) {
    if (!vtk_mode)
        return add_dataset(path, nt, file_dims(sizes_x_first, ncomp), addr, bytes);

    if (path.empty() || path.find('/') != std::string::npos)
        hila::error("hdf5_file " + filename + ": VTKHDF dataset name '" + path +
                    "' must be a plain name without groups");

    // more than 3 directions: drop directions of size 1
    std::vector<uint64_t> sizes;
    std::vector<int64_t> origin;
    for (size_t i = 0; i < sizes_x_first.size(); i++) {
        if (sizes_x_first.size() <= 3 || sizes_x_first[i] > 1) {
            sizes.push_back(sizes_x_first[i]);
            origin.push_back(origin_x_first[i]);
        }
    }
    if (sizes.size() > 3)
        hila::error("hdf5_file " + filename + ": VTKHDF dataset " + path + " has " +
                    std::to_string(sizes.size()) +
                    " dimensions, at most 3 allowed (write a slice, subvolume or projection)");
    while (sizes.size() < 3) {
        sizes.push_back(1);
        origin.push_back(0);
    }

    std::vector<int64_t> extent;
    std::vector<double> orig;
    for (int i = 0; i < 3; i++) {
        extent.push_back(0);
        extent.push_back((int64_t)sizes[i] - 1);
        orig.push_back(origin[i]);
    }

    if (!vtk_grid_set) {
        vtk_grid_set = true;
        vtk_extent = extent;
        vtk_origin = orig;
        object_ref g = group("VTKHDF");
        g.attribute("Version", std::vector<int64_t>{1, 0})
            .attribute("Type", "ImageData")
            .attribute("WholeExtent", vtk_extent)
            .attribute("Origin", vtk_origin)
            .attribute("Spacing", std::vector<double>{1.0, 1.0, 1.0})
            .attribute("Direction", std::vector<double>{1, 0, 0, 0, 1, 0, 0, 0, 1});
    } else if (extent != vtk_extent || orig != vtk_origin) {
        hila::error("hdf5_file " + filename + ": VTKHDF dataset " + path +
                    " has a different grid (size or origin) than the earlier datasets");
    }

    // dataset [nz, ny, nx, ncomp], ncomp omitted if 1
    std::vector<uint64_t> dims = {sizes[2], sizes[1], sizes[0]};
    if (ncomp > 1)
        dims.push_back(ncomp);
    return add_dataset("VTKHDF/PointData/" + path, nt, dims, addr, bytes);
}

/// Write bytes at addr from this rank only (not collective); returns MPI error code
int hdf5_file::write_at(uint64_t addr, const void *buf, uint64_t bytes) {
    int err = MPI_SUCCESS;
    constexpr uint64_t chunk = 1ULL << 30;
    const char *p = static_cast<const char *>(buf);
    for (uint64_t done = 0; done < bytes && err == MPI_SUCCESS; done += chunk) {
        int n = (int)std::min(chunk, bytes - done);
        err = MPI_File_write_at(fh, addr + done, p + done, n, MPI_BYTE, MPI_STATUS_IGNORE);
    }
    return err;
}

/// Write data from rank 0.  Collective
void hdf5_file::write_rank0(uint64_t addr, const void *buf, uint64_t bytes) {
    int err = MPI_SUCCESS;
    if (hila::myrank() == 0)
        err = write_at(addr, buf, bytes);
    hila::io::check_mpi_io(err, "writing", filename);
}

void hdf5_file::add_attribute_raw(int obj, const std::string &name, number_type nt,
                                  const void *data, const std::vector<uint64_t> &dims,
                                  bool scalar) {
    uint64_t n = 1;
    for (auto d : dims)
        n *= d;

    attribute_t a;
    a.name = name;
    a.dtype = encode_datatype(nt);
    a.dspace = encode_dataspace(dims, scalar);
    const uint8_t *p = static_cast<const uint8_t *>(data);
    a.data.assign(p, p + n * nt.size);
    objects[obj].attrs.push_back(a);
}

void hdf5_file::add_attribute(int obj, const std::string &name, const std::string &value) {
    attribute_t a;
    a.name = name;
    // HDF5 strings have size >= 1: the empty string is a single null (padding) byte
    a.dtype = encode_string_datatype(std::max<size_t>(value.size(), 1));
    a.dspace = encode_dataspace(std::vector<uint64_t>(), true);
    a.data.assign(value.begin(), value.end());
    if (value.empty())
        a.data.push_back(0);
    objects[obj].attrs.push_back(a);
}

/// Serialize object headers, children before parents.  header_addr gets the file address
/// of each object header.
void hdf5_file::serialize(int obj, std::vector<uint8_t> &meta, uint64_t meta_start,
                          std::vector<uint64_t> &header_addr) const {
    const object_t &o = objects[obj];
    std::vector<std::pair<uint8_t, bytes_t>> messages;

    if (o.is_group) {
        for (int c : o.children)
            serialize(c, meta, meta_start, header_addr);

        bytes_t linfo;
        put8(linfo, 0); // version
        put8(linfo, 0); // flags: no creation order
        put64(linfo, undefined_address); // fractal heap: none, links stored compactly
        put64(linfo, undefined_address); // name index B-tree: none
        messages.push_back({msg_link_info, linfo});

        bytes_t ginfo;
        put8(ginfo, 0); // version
        put8(ginfo, 0); // flags
        messages.push_back({msg_group_info, ginfo});

        for (int c : o.children) {
            const std::string &name = objects[c].name;
            bytes_t link;
            put8(link, 1); // version
            if (name.size() < 256) {
                put8(link, 0); // flags: hard link, 1-byte name length
                put8(link, name.size());
            } else {
                put8(link, 1); // 2-byte name length
                put16(link, name.size());
            }
            link.insert(link.end(), name.begin(), name.end());
            put64(link, header_addr[c]);
            messages.push_back({msg_link, link});
        }
    } else {
        messages.push_back({msg_dataspace, o.dspace});
        messages.push_back({msg_datatype, o.dtype});

        bytes_t fill;
        put8(fill, 3);    // version
        put8(fill, 0x09); // space allocation early, fill value written if set by user
        messages.push_back({msg_fill_value, fill});

        bytes_t layout;
        put8(layout, 3); // version
        put8(layout, 1); // contiguous
        put64(layout, o.addr);
        put64(layout, o.size);
        messages.push_back({msg_layout, layout});
    }

    for (auto &a : o.attrs)
        messages.push_back({msg_attribute, encode_attribute(a.name, a.dtype, a.dspace, a.data)});

    header_addr[obj] = meta_start + meta.size();
    append(meta, object_header(messages, o.is_group ? "group " + o.name : "dataset " + o.name));
}

void hdf5_file::close() {
    if (!is_open)
        return;

    int err = MPI_SUCCESS;
    if (hila::myrank() == 0) {
        std::vector<uint8_t> meta;
        std::vector<uint64_t> header_addr(objects.size());
        const uint64_t meta_start = (eof + 7) / 8 * 8;
        serialize(0, meta, meta_start, header_addr);

        err = MPI_File_write_at(fh, meta_start, meta.data(), (int)meta.size(), MPI_BYTE,
                                MPI_STATUS_IGNORE);

        bytes_t sb;
        const char signature[] = "\x89HDF\r\n\x1a\n";
        for (int i = 0; i < 8; i++)
            put8(sb, (uint8_t)signature[i]);
        put8(sb, 2);                  // superblock version
        put8(sb, 8);                  // size of offsets
        put8(sb, 8);                  // size of lengths
        put8(sb, 0);                  // file consistency flags
        put64(sb, 0);                 // base address
        put64(sb, undefined_address); // superblock extension: none
        put64(sb, meta_start + meta.size()); // end of file address
        put64(sb, header_addr[0]);    // root group object header
        put32(sb, lookup3(sb.data(), sb.size()));

        if (err == MPI_SUCCESS)
            err = MPI_File_write_at(fh, 0, sb.data(), (int)sb.size(), MPI_BYTE,
                                    MPI_STATUS_IGNORE);
    }
    hila::io::check_mpi_io(err, "writing metadata of", filename);

    hila::io::close_mpi_file(fh, filename);
    is_open = false;
}

} // namespace hila
