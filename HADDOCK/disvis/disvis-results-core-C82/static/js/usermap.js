// Load the Visualization API and the piechart package.
google.charts.load('current', {
    'packages': ['geochart'],
    'mapsApiKey': 'AIzaSyD-9tSrke72PouQMnMX-a7eZSW0jkFMBWY'
});


// Set a callback to run when the Google Visualization API is loaded.
google.charts.setOnLoadCallback(drawMap);


var options = {
    colorAxis: {
        colors: ['#c2d1f0', '#4775d1', '#2952a3', '#142952']
    },
    backgroundColor: '#fff',
    datalessRegionColor: '#ffffff',
    legend: false,
    defaultColor: "#ffffff"
};

var data;
var map;
var table;
var map_data;
var mapdataAPI = $("#mapdataAPI").data().name;


function drawMap() {

    // ------------------------------------------------------------
    // Prepare data for the GeoChat
    let map_data = $.ajax({
        url: mapdataAPI,
        dataType: "json",
        async: false
    }).responseText;
    var obj = JSON.parse(map_data);

    var ready_map_data = Object.keys(obj).map((key) => obj[key]);
    ready_map_data.unshift(["Country", "TotalUsers"]);

    var data = google.visualization.arrayToDataTable(ready_map_data);

    // ------------------------------------------------------------
    // Test data for the GeoChart
    // var data = google.visualization.arrayToDataTable([
    //     ["Country", "TotalUsers"],
    //     ['SY', 500],
    //     ['RU', 600],
    //     ['BR', 700]
    // ]);
    // ------------------------------------------------------------

    map_data = new google.visualization.DataView(data);

    var chart = new google.visualization.GeoChart(document.getElementById('usermap'));
    chart.draw(data, options)
}

// The table is provided by <datatables.net>
$(document).ready(function () {
    $('#table_id').DataTable({
        "lengthMenu": [
            [10, 25, 50, -1],
            [10, 25, 50, "All"]
        ]
    });
});