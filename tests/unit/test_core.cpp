#include <bwa/core/arena.hpp>
#include <bwa/core/vector.hpp>
#include <bwa/core/string.hpp>
#include <bwa/core/hash_map.hpp>
#include <bwa/core/sort.hpp>
#include <iostream>
#include <string>
#include <vector>

using namespace bwa::core;

int main() {
    int passed = 0, failed = 0;

    auto test = [&](const char* name, bool cond) {
        if (cond) { std::cout << "  PASS: " << name << "\n"; ++passed; }
        else { std::cout << "  FAIL: " << name << "\n"; ++failed; }
    };

    // Arena tests
    {
        memory::Arena arena(1024);
        int* arr = arena.allocate<int>(100);
        test("Arena allocate", arr != nullptr);
        for (int i = 0; i < 100; ++i) arr[i] = i;
        test("Arena write/read", arr[99] == 99);
        arena.reset();
        test("Arena reset", arena.used() == 0);
        arena.allocate<int>(10);
        test("Arena after reset", arena.used() == 10 * sizeof(int));
    }

    // Vector tests
    {
        Vector<int> v;
        test("Vector empty", v.empty());
        v.push_back(1);
        v.push_back(2);
        v.push_back(3);
        test("Vector size", v.size() == 3);
        test("Vector access", v[0] == 1 && v[1] == 2 && v[2] == 3);

        Vector<int> v2 = v;
        test("Vector copy", v2 == v);

        Vector<int> v3 = std::move(v);
        test("Vector move", v3.size() == 3 && v.empty());

        v3.clear();
        test("Vector clear", v3.empty());

        v3.reserve(1000);
        test("Vector reserve", v3.capacity() >= 1000);

        v3.shrink_to_fit();
        test("Vector shrink_to_fit", v3.capacity() == 0);

        // SmallVector test
        SmallVector<int> sv;
        for (int i = 0; i < 20; ++i) sv.push_back(i);
        test("SmallVector inline->heap", sv.size() == 20 && sv.capacity() >= 20);

        // Span interop
        std::span<int> sp = v3.span();
        test("Vector span", sp.size() == v3.size());
    }

    // PmrString tests
    {
        PmrString s;
        test("String empty", s.empty());
        s.kputs("Hello");
        test("String kputs", s.view() == "Hello");
        s.kputc(' ');
        test("String kputc", s.view() == "Hello ");
        s.kputs("World");
        test("String append", s.view() == "Hello World");
        s.kputw(42);
        test("String kputw", s.view() == "Hello World42");
        s.kputuw(100u);
        test("String kputuw", s.view() == "Hello World42100");
        s.kputl(-1234567890L);
        test("String kputl", s.view().find("-1234567890") != std::string_view::npos);
        s.kputul(1234567890ul);
        test("String kputul", s.view().find("1234567890") != std::string_view::npos);
        s.kputd(3.14159);
        test("String kputd", s.view().find("3.14") != std::string_view::npos);
        s.kputf(2.5, 2);
        test("String kputf", s.view().find("2.50") != std::string_view::npos);

        // Format
        s.clear();
        s.ksprintf("Value: {}, hex: {:x}", 255, 255);
        test("String ksprintf", s.view() == "Value: 255, hex: ff");

        // Copy/move
        PmrString s2 = s;
        test("String copy", s2 == s);
        PmrString s3 = std::move(s2);
        test("String move", s3 == s && s2.empty());

        // Span
        std::string_view sv = s3.view();
        test("String string_view", sv == s3.view());

        // Find
        test("String find", s3.find("255") != PmrString::NPOS);
        test("String rfind", s3.rfind("f") != PmrString::NPOS);

        // Substr
        PmrString sub = s3.substr(0, 5);
        test("String substr", sub.view() == "Value");
    }

    // HashMap tests
    {
        PmrHashMap<int, int> m;
        test("HashMap empty", m.empty());
        m.insert(1, 10);
        m.insert(2, 20);
        m.insert(3, 30);
        test("HashMap size", m.size() == 3);
        test("HashMap find", m.find(2)->second == 20);
        test("HashMap contains", m.contains(3) && !m.contains(4));
        test("HashMap operator[]", m[4] == 0);
        m[4] = 40;
        test("HashMap operator[] assign", m[4] == 40);
        test("HashMap at", m.at(1) == 10);
        m.erase(2);
        test("HashMap erase", m.size() == 3 && !m.contains(2));
        m.clear();
        test("HashMap clear", m.empty());

        // String keys
        PmrHashMap<PmrString, int> sm;
        sm.insert("foo", 1);
        sm.insert("bar", 2);
        test("HashMap string key", sm.find("foo")->second == 1);

        // khash compatibility
        sm.kh_put("baz");
        test("HashMap kh_put", sm.contains("baz"));
        sm.kh_del("bar");
        test("HashMap kh_del", !sm.contains("bar"));
    }

    // HashSet tests
    {
        PmrHashSet<int> s;
        s.insert(1); s.insert(2); s.insert(3);
        test("HashSet size", s.size() == 3);
        test("HashSet contains", s.contains(2));
        s.erase(2);
        test("HashSet erase", !s.contains(2));
    }

    // Sort tests
    {
        Vector<int> v;
        for (int i = 100; i >= 0; --i) v.push_back(i);
        core::radix_sort(v);
        bool ok = true;
        for (int i = 0; i <= 100; ++i) if (v[i] != i) ok = false;
        test("RadixSort unsigned", ok);

        Vector<int> v2;
        for (int i = -50; i <= 50; ++i) v2.push_back(i);
        core::radix_sort(v2);
        ok = true;
        for (int i = 0; i <= 100; ++i) if (v2[i] != i - 50) ok = false;
        test("RadixSort signed", ok);

        Vector<float> v3;
        for (int i = 100; i >= 0; --i) v3.push_back(float(i) * 0.1f);
        core::radix_sort(v3);
        ok = true;
        for (int i = 0; i <= 100; ++i) if (v3[i] != float(i) * 0.1f) ok = false;
        test("RadixSort float", ok);

        // Pair sort
        Vector<uint32_t> keys = {3, 1, 4, 1, 5, 9};
        Vector<int> vals = {30, 10, 40, 11, 50, 90};
        core::radix_sort_pairs(keys, vals);
        ok = true;
        for (size_t i = 1; i < keys.size(); ++i) if (keys[i-1] > keys[i]) ok = false;
        test("RadixSort pairs", ok && vals[0] == 10 && vals[1] == 11);

        // Argsort
        Vector<int> v4 = {5, 2, 8, 1, 9};
        auto idx = core::argsort(v4);
        ok = idx.size() == 5 && idx[0] == 3 && idx[1] == 1 && idx[4] == 4;
        test("Argsort", ok);
    }

    // PackedSequence tests
    {
        using index::PackedSequence;
        PackedSequence seq;
        seq.append("ACGTACGTNN", 10);
        test("PackedSequence size", seq.size() == 10);
        test("PackedSequence get", seq.get(0) == 0 && seq.get(1) == 1 &&
             seq.get(2) == 2 && seq.get(3) == 3);
        test("PackedSequence N", seq.get(8) == 4 && seq.is_n(8));
        test("PackedSequence not N", !seq.is_n(0));

        // Iteration
        int count = 0;
        for (uint8_t b : seq) { (void)b; ++count; }
        test("PackedSequence iter", count == 10);

        // Reverse complement
        PackedSequence seq2;
        seq2.append("ACGT", 4);
        seq2.reverse_complement();
        test("PackedSequence revcomp", seq2.get(0) == 3 && seq2.get(1) == 2 &&
             seq2.get(2) == 1 && seq2.get(3) == 0);

        // To string
        PmrString out;
        seq2.to_string(out);
        test("PackedSequence to_string", out.view() == "ACGT");
    }

    std::cout << "\nResults: " << passed << " passed, " << failed << " failed\n";
    return failed ? 1 : 0;
}