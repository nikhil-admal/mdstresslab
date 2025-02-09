
#include <iostream>
#include <string>
#include <initializer_list>
#include <curl/curl.h>

// A callback that writes data into a std::string
static size_t write_callback(void* contents, size_t size, size_t nmemb, void* userp)
{
    const size_t totalSize = size * nmemb;
    std::string& buffer = *static_cast<std::string*>(userp);
    buffer.append(static_cast<char*>(contents), totalSize);
    return totalSize;
}

/**
 * kim_query() with default parameters:
 *   1) query_function        - required
 *   2) model_name            - default is "" (no model)
 *   3) extra_args (init list) - default is empty
 */
std::string kim_query(
    const std::string& query_function,
    const std::initializer_list<std::string>& extra_args ={
      "species=[]",
    }
)
{
    // URL to contact
    const std::string url = "https://query.openkim.org/api/" + query_function;

    // Build the POST fields
    std::string post_data;

    // Append extra_args, each one is "key=[stuff]"
    for (auto& kv : extra_args) {
        // If there's already something in post_data, prepend '&'
        if (!post_data.empty()) {
            post_data += "&";
        }
        post_data += kv;
    }

    // Initialize response string
    std::string response;

    // Initialize libcurl
    curl_global_init(CURL_GLOBAL_DEFAULT);
    CURL* curl_handle = curl_easy_init();
    if (!curl_handle) {
        curl_global_cleanup();
        return "Error: Failed to initialize libcurl!";
    }

    curl_easy_setopt(curl_handle, CURLOPT_SSL_VERIFYHOST, 0);
    curl_easy_setopt(curl_handle, CURLOPT_SSL_VERIFYPEER, 0);

    // Set the options
    curl_easy_setopt(curl_handle, CURLOPT_URL, url.c_str());
    curl_easy_setopt(curl_handle, CURLOPT_POST, 1L);
    curl_easy_setopt(curl_handle, CURLOPT_POSTFIELDS, post_data.c_str());
    curl_easy_setopt(curl_handle, CURLOPT_WRITEFUNCTION, write_callback);
    curl_easy_setopt(curl_handle, CURLOPT_WRITEDATA, &response);

    // Perform the request
    CURLcode res = curl_easy_perform(curl_handle);
    if (res != CURLE_OK) {
        response = std::string("curl error: ") + curl_easy_strerror(res);
    }

    // Cleanup
    curl_easy_cleanup(curl_handle);
    curl_global_cleanup();

    return response;
}

/*
int main()
{

    std::string result2 = kim_query(
        "get_available_models",
        {
          "species=[\"Si\",\"C\"]",
          "species_logic=[\"and\"]",
          "model_interface=[\"sm\"]",
          "potential_type=[\"meam\"]",
          "simulator_name=[\"LAMMPS\"]"
        }
    );
    std::cout << "\nResponse from custom query:\n" << result2 << std::endl;
     std::string result = kim_query(
        "get_lattice_constant_cubic",
        {
            "model=[\"MO_123629422045_005\"]",
            "crystal=[\"fcc\"]",
            "species=[\"Al\"]",
            "units=[\"angstrom\"]"
        }
    );

    std::cout << result << "\n";
    return 0;
}
*/
