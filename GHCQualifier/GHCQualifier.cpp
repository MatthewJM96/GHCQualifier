#include <fstream>
#include <sstream>
#include <string>
#include <iostream>
#include <vector>
#include <algorithm>
#include <random>
#include <limits>
#include <map>
#include <unordered_map>
#include <functional>

#pragma region DataStructures
using ID    = size_t;
using Count = size_t;

using CachedVideos = std::vector<ID>;
using UsedCache    = std::pair<ID, CachedVideos>;
using UsedCaches   = std::unordered_map<ID, CachedVideos>;

using VideoTimes    = std::unordered_map<ID, double>;
using EndpointTimes = std::unordered_map<ID, VideoTimes>;

using CachePriorityList = std::map<double, ID, std::greater<double>>;

using SolverResponse = std::pair<ID, bool>;
using HPCSResponse   = SolverResponse;
using HPVSResponse   = SolverResponse;

struct VideoDatum {
    int   size;
    Count numReq;
};
using VideoData = std::map<ID, VideoDatum>;
using VideoPriorityList = std::map<double, ID, std::greater<double>>;

struct Metadata {
    Count videoCount;
    Count endpointCount;
    Count requestCount;
    Count cacheCount;
    int   cacheSize;
    int   smallestVideoSize;
};

struct Video {
    ID  id;
    int size;
};
using Videos = std::vector<Video>;

struct CacheConnection {
    ID  cacheId;
    int latency;
};
using CacheConnections = std::vector<CacheConnection>;

struct Request {
    ID    id;
    ID    videoId;
    Count requestCount;
};
using Requests = std::vector<Request>;

struct Endpoint {
    ID               id;
    int              latency;
    Count            cacheCount;
    CacheConnections cacheConnections;
    Requests         requests;
};
using Endpoints = std::vector<Endpoint>;

struct Cache {
    ID        id;
    Endpoints connections;
    int       remainingSpace;
};
using Caches = std::vector<Cache>;

struct Data {
    Metadata  metadata;
    Videos    videos;
    Endpoints endpoints;
    Caches    caches;
};
#pragma endregion

/// Model that acquires the data.
class Loader {
public:
    Loader() : m_data({}) {}
    ~Loader() {
        dispose();
    }

    void init(std::string filepath = "example.in") {
        m_filepath = filepath;
    }
    void dispose() {
        // Close file is still open.
        if (m_file.is_open()) {
            m_file.close();
        }

        // Clear up data.
        m_data = {};
    }

    void setFilepath(std::string filepath) {
        m_filepath = filepath;
    }

    Data getData() {
        // If we haven't yet loaded data from the file, do so.
        if (m_data.metadata.videoCount == 0) {
            loadDataFromFile();
        }
        return m_data;
    }
private:
    bool openFile() {
        m_file.open(m_filepath, std::ios::in);
        return m_file.is_open();
    }

    void loadDataFromFile() {
        // If file isn't already open, and failed to open on an attempt, exit the program.
        if (!m_file.is_open() && !openFile()) {
            std::cout << "Could not open file: " << m_filepath << "." << std::endl
                << "Exiting..." << std::endl;
            std::getchar();
            exit(0);
        }

        std::string line;

        // Grab metadata.
        std::getline(m_file, line);
        std::stringstream mdStream(line);

        mdStream >> m_data.metadata.videoCount;
        mdStream >> m_data.metadata.endpointCount;
        mdStream >> m_data.metadata.requestCount;
        mdStream >> m_data.metadata.cacheCount;
        mdStream >> m_data.metadata.cacheSize;

        // Grab videos.
        std::getline(m_file, line);
        std::stringstream vStream(line);

        int smallestVideoSize = INT_MAX;
        for (size_t i = 0; i < m_data.metadata.videoCount; ++i) {
            int videoSize;
            vStream >> videoSize;

            if (smallestVideoSize > videoSize) {
                smallestVideoSize = videoSize;
            }

            m_data.videos.push_back({ i, videoSize });
        }
        m_data.metadata.smallestVideoSize = smallestVideoSize;

        // Grab endpoints.
        for (size_t i = 0; i < m_data.metadata.endpointCount; ++i) {
            std::getline(m_file, line);
            std::stringstream eStream(line);

            int datacenterLatency;
            size_t cacheConnectionsCount;
            CacheConnections cacheConnections;

            eStream >> datacenterLatency;
            eStream >> cacheConnectionsCount;

            for (size_t j = 0; j < cacheConnectionsCount; ++j) {
                std::getline(m_file, line);
                std::stringstream ecStream(line);

                size_t cacheId;
                int cacheLatency;
                ecStream >> cacheId;
                ecStream >> cacheLatency;

                cacheConnections.push_back({ cacheId, cacheLatency });
            }

            m_data.endpoints.push_back({ i, datacenterLatency, cacheConnectionsCount, cacheConnections });
        }

        // Grab requests.
        for (size_t i = 0; i < m_data.metadata.requestCount; ++i) {
            std::getline(m_file, line);
            std::stringstream rStream(line);

            size_t videoId, endpointId, requestCount;
            rStream >> videoId;
            rStream >> endpointId;
            rStream >> requestCount;

            m_data.endpoints[endpointId].requests.push_back({ i, videoId, requestCount });
        }

        constructCaches();
    }

    // Constructs list of data from cache point of view.
    void constructCaches() {
        m_data.caches.resize(m_data.metadata.cacheCount);
        for (size_t i = 0; i < m_data.metadata.cacheCount; ++i) {
            m_data.caches[i] = { i, {}, m_data.metadata.cacheSize };
        }
        
        for (Endpoint endpoint : m_data.endpoints) {
            for (CacheConnection cacheConnection : endpoint.cacheConnections) {
                m_data.caches[cacheConnection.cacheId].connections.push_back({ endpoint.id, endpoint.latency - cacheConnection.latency, 0, {}, endpoint.requests });
            }
        }
    }

    std::fstream m_file;
    std::string m_filepath;

    Data m_data;
};

/// Interface for solvers.
template <typename Result>
class ISolver {
public:
    virtual void solve() = 0;

    Result getResult() const {
        return m_result;
    }
protected:
    Result m_result;
};

// Generic simulated annealing implementation.
template <typename State, typename Energy = double, typename Temperature = double,
    typename Generator = std::mt19937_64, typename Count = size_t,
    typename std::enable_if<std::is_arithmetic<Energy>::value && std::is_arithmetic<Temperature>::value && std::is_arithmetic<Count>::value>::type* = 0>
    class SimulatedAnnealer {
    public:
        using EnergyFunc = Energy(*)(State);
        using TemperatureFunc = Temperature(*)(Count);
        using NextFunc = State(*)(State, Energy);

        SimulatedAnnealer() {}
        ~SimulatedAnnealer() {
            dispose();
        }

        SimulatedAnnealer* init(EnergyFunc eFunc = nullptr, TemperatureFunc tFunc = nullptr, NextFunc nFunc = nullptr) {
            if (m_initialised) return this;
            m_initialised = true;

            m_energyFunc = eFunc;
            m_temperatureFunc = tFunc;
            m_nextFunc = nFunc;

            return this;
        }
        SimulatedAnnealer* dispose() {
            if (!m_initialised) return this;
            m_initialised = false;

            m_energyFunc = nullptr;
            m_temperatureFunc = nullptr;
            m_nextFunc = nullptr;

            return this;
        }

        SimulatedAnnealer* setEnergyFunc(EnergyFunc eFunc) {
            m_energyFunc = eFunc;
            return this;
        }
        SimulatedAnnealer* setTemperatureFunc(TemperatureFunc tFunc) {
            m_temperatureFunc = tFunc;
            return this;
        }
        SimulatedAnnealer* setNextFunc(NextFunc nFunc) {
            m_nextFunc = nFunc;
            return this;
        }

        State run(State initialState, Count iterationCount) {
            std::random_device rd;
            Generator generator(rd());
#define oldState initialState
#define count iterationCount
            Energy oldEnergy = m_energyFunc(initialState);

            State  bestState = initialState;
            Energy bestEnergy = oldEnergy;

            std::uniform_real_distribution<Energy> randFunc((Energy)0, (Energy)1);

            for (; count > 0; --count) {
                State newState = m_nextFunc(oldState, std::abs(bestEnergy - oldEnergy));
                Energy newEnergy = m_energyFunc(newState);

                if (newEnergy < bestEnergy) {
                    bestState = newState;
                    bestEnergy = newEnergy;
                    oldState = std::move(newState);
                    oldEnergy = std::move(newEnergy);
                    continue;
                }

                Temperature t = m_temperatureFunc(count);
                if (newEnergy < oldEnergy || std::exp((oldEnergy - newEnergy) / t) > randFunc(generator)) {
                    oldState = std::move(newState);
                    oldEnergy = std::move(newEnergy);
                }
            }

            return bestState;
#undef oldState
#undef count
        }
    private:
        bool m_initialised;

        EnergyFunc      m_energyFunc;
        TemperatureFunc m_temperatureFunc;
        NextFunc        m_nextFunc;
};

EndpointTimes getTimeToNearestStore(Data data, UsedCaches usedCaches) {
    EndpointTimes et;
    for (Endpoint endpoint : data.endpoints) {
        // Insert endpoint into nearest times var.
        et.insert({ endpoint.id, {} });
        // For each video request at this endpoint, add it to the nearest times var with time to nearest store of the requested video.
        for (Request request : endpoint.requests) {
            et[endpoint.id].insert({ request.videoId, 0.0 });

            double bestTime = endpoint.latency; // Set best time by default to be to the datacenter.
            for (auto usedCache : usedCaches) {
                // If the video exists in any of the used caches, then see if that cache is connected to this endpoint.
                if (std::find(usedCache.second.begin(), usedCache.second.end(), request.videoId) != usedCache.second.end()) {
                    auto& it = std::find_if(endpoint.cacheConnections.begin(), endpoint.cacheConnections.end(), [&usedCache](CacheConnection a) {
                        return a.cacheId == usedCache.first;
                    });
                    // If a cache storing the video is connected to this endpoint, and its latency is better than the current best latency to fetch the video, set the latency of this endpoint -> this cache as the new best latency.
                    if (it != endpoint.cacheConnections.end()) {
                        if (bestTime > it->latency) {
                            bestTime = it->latency;
                        }
                    }
                }
            }

            // Insert found best latency into nearest times var.
            et[endpoint.id][request.videoId] = bestTime;
        }
    }
    return et;
}

class HighestPriorityCacheSolver : public ISolver<HPCSResponse> {
public:
    ~HighestPriorityCacheSolver() {
        dispose();
    }

    void init(Data data, EndpointTimes nearestTimes, UsedCaches usedCaches, double weighting) {
        m_k            = weighting;
        m_nearestTimes = nearestTimes;
        m_usedCaches   = usedCaches;
        m_data         = data;
        m_skipCount    = 0;
    }
    void dispose() {
        EndpointTimes().swap(m_nearestTimes);
        UsedCaches().swap(m_usedCaches);
    }

    void setSkipCount(Count skips) {
        m_skipCount = skips;
    }

    void solve() {
        CachePriorityList cpl;

        for (Cache cache : m_data.caches) {
            cpl.insert({ calculateGoodness(cache), cache.id });
        }

        m_result = HPCSResponse(0, false);
        Count skips = 0;
        auto& it = cpl.begin();
        while (it != cpl.end()) {
            if (m_data.caches[it->second].remainingSpace > m_data.metadata.smallestVideoSize) {
                if (skips != m_skipCount) {
                    ++it;
                    ++skips;
                    continue;
                }
                m_result = HPCSResponse(it->second, true);
                break;
            }
            ++it;
        }
    }
private:
    double calculateGoodness(Cache cache) {
        double res = (double)cache.remainingSpace;

        for (Endpoint connection : cache.connections) {
            for (Request request : connection.requests) {
                // If the cache we are considering already has the video, skip.
                auto& usedCache = std::find_if(m_usedCaches.begin(), m_usedCaches.end(), [&cache](const UsedCache& usedCache) {
                    return usedCache.first == cache.id;
                });
                if (usedCache != m_usedCaches.end() && std::find(usedCache->second.begin(), usedCache->second.end(), request.videoId) != usedCache->second.end()) continue;

                // Calculate goodness.
                res += m_k * m_nearestTimes[connection.id][request.videoId] * (double)request.requestCount / (double)m_data.videos[request.videoId].size;
            }
        }

        return res;
    }

    Data          m_data;
    EndpointTimes m_nearestTimes;
    UsedCaches    m_usedCaches;
    double        m_k;
    Count         m_skipCount;
};

class HighestPriorityVideoSolver : public ISolver<HPVSResponse> {
    // Constructs priority lists for videos of a given cache.
public:
    ~HighestPriorityVideoSolver() {
        dispose();
    }

    void init(Data data, EndpointTimes nearestTimes, UsedCaches usedCaches, ID id) {
        m_nearestTimes = nearestTimes;
        m_usedCaches = usedCaches;
        m_data = data;
        m_id = id;
        m_skipCount = 0;
    }
    void dispose() {
        EndpointTimes().swap(m_nearestTimes);
        UsedCaches().swap(m_usedCaches);
    }

    void setSkipCount(Count skips) {
        m_skipCount = skips;
    }

    void solve() {
        // Construct video data object, we use this to calculate goodness later.
        VideoData vd;
        for (Endpoint connection : m_data.caches[m_id].connections) {
            for (Request request : connection.requests) {
                // If the cache we are considering already has the video, skip.
                auto& usedCache = std::find_if(m_usedCaches.begin(), m_usedCaches.end(), [&](const UsedCache& usedCache) {
                    return usedCache.first == m_id;
                });
                if (usedCache != m_usedCaches.end() && std::find(usedCache->second.begin(), usedCache->second.end(), request.videoId) != usedCache->second.end()) continue;

                auto& it = vd.find(request.videoId);
                if (it == vd.end()) {
                    vd.insert({ request.videoId, { m_data.videos[request.videoId].size, 0 } });
                }
                vd.at(request.videoId).numReq += request.requestCount;
            }
        }

        // Kinda sucks but ceebs fixing for now.
        std::unordered_map<ID, double> temp;
        for (Endpoint connection : m_data.caches[m_id].connections) {
            for (Request request : connection.requests) {
                auto& datumIt = vd.find(request.videoId);
                if (datumIt == vd.end()) continue;

                auto& it = temp.find(request.videoId);
                if (it == temp.end()) {
                    temp.insert({ request.videoId, 0.0 });
                }

                temp.at(request.videoId) += m_nearestTimes[connection.id][request.videoId] * (double)datumIt->second.numReq / (double)datumIt->second.size;
            }
        }

        VideoPriorityList vpl;
        for (auto& it = temp.begin(); it != temp.end(); ++it) {
            vpl.insert({ it->second, it->first });
        }

        m_result = HPVSResponse(0, false);
        Count skips = 0;
        auto& it = vpl.begin();
        while (it != vpl.end()) {
            if (m_data.caches[m_id].remainingSpace > m_data.videos[it->second].size) {
                if (skips != m_skipCount) {
                    ++it;
                    ++skips;
                    continue;
                }
                m_result = HPVSResponse(it->second, true);
                break;
            }
            ++it;
        }
    }
private:
    Data          m_data;
    EndpointTimes m_nearestTimes;
    UsedCaches    m_usedCaches;
    size_t        m_id;
    Count         m_skipCount;
};

bool cachesFull(Data data) {
    for (auto& cache : data.caches) {
        if (cache.remainingSpace > data.metadata.smallestVideoSize) return false;
    }
    return true;
}

void doThatShit(std::string filename) {
    Loader loader;
    loader.init(filename + ".in");

    Data data = loader.getData();

    UsedCaches usedCaches;

    do {
        EndpointTimes nt = getTimeToNearestStore(data, usedCaches);

        HighestPriorityCacheSolver hpcs;
        hpcs.init(data, nt, usedCaches, 0.0001);


        HPCSResponse cacheResponse;
        HPVSResponse videoResponse;
        bool  skipMaxExceededForCache = false;
        bool  skipMaxExceededForVideo = false;
        Count skipsForCache = 0;
        Count skipsForVideo = 0;
        do {
            skipMaxExceededForVideo = false;

            hpcs.setSkipCount(skipsForCache);
            hpcs.solve();

            cacheResponse = hpcs.getResult();
            if (!cacheResponse.second) {
                ++skipsForCache;
                if (skipsForCache > 5) {
                    skipMaxExceededForCache = true;
                    break;
                }
                continue;
            }

            HighestPriorityVideoSolver hpvs;
            hpvs.init(data, nt, usedCaches, cacheResponse.first);

            while (true) {
                hpvs.setSkipCount(skipsForVideo);
                hpvs.solve();

                videoResponse = hpvs.getResult();
                if (videoResponse.second) break;

                ++skipsForVideo;
                if (skipsForVideo > 5) {
                    skipMaxExceededForVideo = true;
                    break;
                }
            }
            hpvs.dispose();

            if (skipMaxExceededForVideo) {
                ++skipsForCache;
                if (skipsForCache > 5) {
                    skipMaxExceededForCache = true;
                    break;
                }
            }
        } while (!cacheResponse.second || !videoResponse.second);
        hpcs.dispose();

        if (skipMaxExceededForCache || !videoResponse.second) break;
        //if (skipMaxExceededForCache) break;

        if (usedCaches.find(cacheResponse.first) == usedCaches.end()) {
            usedCaches.insert({ cacheResponse.first,{} });
        }
        usedCaches[cacheResponse.first].push_back(videoResponse.first);

        data.caches[cacheResponse.first].remainingSpace -= data.videos[videoResponse.first].size;
    } while (!cachesFull(data));

    std::fstream file(filename + ".txt", std::ios::out);
    file << usedCaches.size() << "\n";
    for (auto usedCache : usedCaches) {
        file << usedCache.first;
        for (auto i : usedCache.second) {
            file << " " << i;
        }
        file << "\n";
    }
    file.close();
}

int main() {
    doThatShit("me_at_the_zoo");
    doThatShit("videos_worth_spreading");

    std::cout << "Press any key to exit..." << std::endl;
    getchar();
    return 0;
}