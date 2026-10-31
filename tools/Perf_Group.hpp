
#ifndef KACHE_HASH_PERF_GROUP_HPP
#define KACHE_HASH_PERF_GROUP_HPP



#include <cstdint>
#include <cstddef>
#include <string>
#include <vector>
#include <cstring>
#include <iostream>
#include <cerrno>
#include <cstdlib>
#include <algorithm>

#ifdef __linux__
#include <linux/perf_event.h>
#include <sys/ioctl.h>
#include <sys/syscall.h>
#include <unistd.h>
#endif


// A group of performance-counter events, counting the user-space events of the
// constructing thread and of the threads spawned afterwards; hence it is to be
// constructed before spawning any thread. Cycles, instructions, and a few
// software events are always counted alongside the requested raw hardware
// events. Counting is supported on Linux only.
class Perf_Group
{
private:

    // An event to count.
    struct Event
    {
        std::string name;   // Name of the event.
        uint32_t type;  // Type of the event, per `perf_event_attr`.
        uint64_t config;    // Configuration of the event, per `perf_event_attr`.
        bool grouped;   // Whether the event is scheduled together with the hardware group.
        int fd; // File descriptor of the event.
        uint64_t val[3];    // Count, time enabled, and time running as of the last `start` or `stop`.
        uint64_t delta[3];  // Count, time enabled, and time running between the last `start` and `stop`.

        Event(const std::string& name, const uint32_t type, const uint64_t config, const bool grouped):
              name(name), type(type), config(config), grouped(grouped), fd(-1), val(), delta()
        {}
    };

    std::vector<Event> E;   // The events; `E[0]` leads the hardware group.


    // Puts the count, time enabled, and time running of the event `e` so far at
    // `r`.
    static void read_event(const Event& e, uint64_t (&r)[3]);

public:

    // Constructs a group of the raw hardware events in `spec`, a comma-separated
    // list of `name=rCODE` items with hexadecimal raw codes `CODE`, alongside the
    // fixed events. Exits on failure.
    Perf_Group(const std::string& spec);

    Perf_Group(const Perf_Group&) = delete;
    Perf_Group& operator=(const Perf_Group&) = delete;

    ~Perf_Group();

    // Starts counting.
    void start();

    // Stops counting.
    void stop();

    // Returns the number of events.
    std::size_t size() const { return E.size(); }

    // Returns the name of the `i`'th event.
    const std::string& name(const std::size_t i) const { return E[i].name; }

    // Returns the count of the `i`'th event between the last `start` and `stop`.
    uint64_t count(const std::size_t i) const { return E[i].delta[0]; }

    // Returns the percentage of the time between the last `start` and `stop` that
    // the hardware group was counting for; it is less than 100 iff the group was
    // multiplexed with other events.
    double running_pct() const { return E[0].delta[1] > 0 ? 100.0 * E[0].delta[2] / E[0].delta[1] : 0; }
};


inline Perf_Group::Perf_Group(const std::string& spec)
{
#ifdef __linux__
    E.emplace_back("cycles", PERF_TYPE_HARDWARE, PERF_COUNT_HW_CPU_CYCLES, true);
    E.emplace_back("instructions", PERF_TYPE_HARDWARE, PERF_COUNT_HW_INSTRUCTIONS, true);

    for(std::size_t b = 0; b < spec.size(); )
    {
        const auto e = std::min(spec.find(',', b), spec.size());
        const auto item = spec.substr(b, e - b);
        const auto eq = item.find('=');
        if(eq == std::string::npos || eq + 2 >= item.size() || item[eq + 1] != 'r')
        {
            std::cerr << "Invalid raw event " << item << "; expected `name=rCODE`. Aborting.\n";
            std::exit(EXIT_FAILURE);
        }

        E.emplace_back(item.substr(0, eq), PERF_TYPE_RAW, std::strtoull(item.c_str() + eq + 2, nullptr, 16), true);
        b = e + 1;
    }

    E.emplace_back("task_clock_ns", PERF_TYPE_SOFTWARE, PERF_COUNT_SW_TASK_CLOCK, false);
    E.emplace_back("page_faults", PERF_TYPE_SOFTWARE, PERF_COUNT_SW_PAGE_FAULTS, false);
    E.emplace_back("context_switches", PERF_TYPE_SOFTWARE, PERF_COUNT_SW_CONTEXT_SWITCHES, false);
    E.emplace_back("cpu_migrations", PERF_TYPE_SOFTWARE, PERF_COUNT_SW_CPU_MIGRATIONS, false);

    for(auto& e : E)
    {
        const bool leader = (!e.grouped || &e == &E[0]);

        perf_event_attr attr;
        std::memset(&attr, 0, sizeof(attr));
        attr.size = sizeof(attr);
        attr.type = e.type;
        attr.config = e.config;
        attr.disabled = leader; // Members of the group count whenever their leader does.
        attr.inherit = 1;   // Counts the threads spawned afterwards too.
        attr.exclude_kernel = 1;
        attr.exclude_hv = 1;
        attr.read_format = PERF_FORMAT_TOTAL_TIME_ENABLED | PERF_FORMAT_TOTAL_TIME_RUNNING;

        e.fd = syscall(SYS_perf_event_open, &attr, 0, -1, (leader ? -1 : E[0].fd), PERF_FLAG_FD_CLOEXEC);
        if(e.fd < 0)
        {
            std::cerr << "Error opening the performance-counter event " << e.name << ": " << std::strerror(errno) << ". Aborting.\n";
            std::exit(EXIT_FAILURE);
        }
    }
#else
    (void)spec;
    std::cerr << "Performance counters are supported on Linux only. Aborting.\n";
    std::exit(EXIT_FAILURE);
#endif
}


inline Perf_Group::~Perf_Group()
{
#ifdef __linux__
    for(auto it = E.rbegin(); it != E.rend(); ++it)
        close(it->fd);
#endif
}


inline void Perf_Group::read_event(const Event& e, uint64_t (&r)[3])
{
#ifdef __linux__
    if(read(e.fd, r, sizeof(r)) != sizeof(r))
    {
        std::cerr << "Error reading the performance-counter event " << e.name << ": " << std::strerror(errno) << ". Aborting.\n";
        std::exit(EXIT_FAILURE);
    }
#else
    (void)e, (void)r;
#endif
}


inline void Perf_Group::start()
{
#ifdef __linux__
    for(auto& e : E)
        read_event(e, e.val);   // Counts from exited threads persist, so counts are taken as differences.

    ioctl(E[0].fd, PERF_EVENT_IOC_ENABLE, PERF_IOC_FLAG_GROUP);
    for(const auto& e : E)
        if(!e.grouped)
            ioctl(e.fd, PERF_EVENT_IOC_ENABLE, 0);
#endif
}


inline void Perf_Group::stop()
{
#ifdef __linux__
    ioctl(E[0].fd, PERF_EVENT_IOC_DISABLE, PERF_IOC_FLAG_GROUP);
    for(const auto& e : E)
        if(!e.grouped)
            ioctl(e.fd, PERF_EVENT_IOC_DISABLE, 0);

    for(auto& e : E)
    {
        uint64_t r[3];
        read_event(e, r);
        for(std::size_t i = 0; i < 3; ++i)
            e.delta[i] = r[i] - e.val[i], e.val[i] = r[i];
    }
#endif
}



#endif
