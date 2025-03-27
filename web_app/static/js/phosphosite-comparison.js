/**
 * Modified Phosphosite Comparison Visualization
 * This script creates interactive visualizations comparing known vs unknown phosphosites.
 * Auto-initializes without requiring toggle button clicks.
 */

// Configuration for visualization metrics and colors remains the same
const COMPARISON_CONFIG = {
    colors: {
      known: '#1976d2',      // Bright blue for known sites
      unknown: '#ff9800',    // Orange for unknown sites
      gradient: ['#e3f2fd', '#bbdefb', '#90caf9', '#64b5f6', '#42a5f5', '#2196f3', '#1e88e5', '#1976d2', '#1565c0', '#0d47a1'],
      categorical: ['#1976d2', '#ff9800', '#4caf50', '#f44336', '#9c27b0', '#ff5722', '#607d8b']
    },
    metrics: [
      {
        id: 'plddt',
        name: 'Mean pLDDT Score',
        description: 'Higher pLDDT scores indicate greater confidence in the protein structure prediction.',
        format: (value) => value.toFixed(1)
      },
      {
        id: 'site_plddt',
        name: 'Site pLDDT Score',
        description: 'pLDDT score specifically at the phosphorylation site position.',
        format: (value) => value.toFixed(1)
      },
      {
        id: 'nearby',
        name: 'Nearby Residues (10Å)',
        description: 'Number of residues within 10Å of the phosphorylation site. Higher values suggest a more buried position.',
        format: (value) => value.toFixed(0)
      },
      {
        id: 'surface',
        name: 'Surface Accessibility (%)',
        description: 'Percentage of the site that is accessible to solvent. Higher values indicate greater exposure.',
        format: (value) => value.toFixed(1) + '%'
      },
      {
        id: 'acidic',
        name: 'Acidic Residue Content (%)',
        description: 'Percentage of acidic residues (D, E) in the vicinity of the phosphosite.',
        format: (value) => value.toFixed(1) + '%'
      },
      {
        id: 'basic',
        name: 'Basic Residue Content (%)',
        description: 'Percentage of basic residues (K, R, H) in the vicinity of the phosphosite.',
        format: (value) => value.toFixed(1) + '%'
      }
    ]
  };
  
  // Modified initialization function
  function initPhosphositeComparison(containerId) {
    console.log(`Initializing phosphosite comparison in container: ${containerId}`);
    const container = document.getElementById(containerId);
    
    if (!container) {
      console.error(`Container element with ID "${containerId}" not found.`);
      return;
    }
    
    // Extract phosphosite data from the table
    const phosphosites = extractPhosphositesFromTable();
    console.log(`Extracted ${phosphosites.length} phosphosites for comparison`);
    
    if (!phosphosites || phosphosites.length === 0) {
      container.innerHTML = '<div class="alert alert-warning">No phosphosite data available for comparison analysis.</div>';
      return;
    }
    
    // Split the phosphosites into known and unknown groups
    const knownSites = phosphosites.filter(site => site.isKnown);
    const unknownSites = phosphosites.filter(site => !site.isKnown);
    
    console.log(`Found ${knownSites.length} known phosphosites and ${unknownSites.length} unknown phosphosites`);
    
    if (knownSites.length === 0 || unknownSites.length === 0) {
      container.innerHTML = '<div class="alert alert-warning">Unable to perform comparison: need both known and unknown phosphosites.</div>';
      return;
    }
    
    // Create the container structure for the visualizations
    container.innerHTML = `
      <div class="card">
        <div class="card-header">
          <h5 class="mb-0">Phosphosite Profile Comparison: Known vs Unknown Sites</h5>
        </div>
        <div class="card-body">
          <div class="row mb-4">
            <div class="col-12">
              <div class="d-flex align-items-center justify-content-center mb-3">
                <div class="d-flex align-items-center me-4">
                  <div style="width: 16px; height: 16px; background-color: ${COMPARISON_CONFIG.colors.known}; border-radius: 3px; margin-right: 6px;"></div>
                  <span>Known Phosphosites (${knownSites.length})</span>
                </div>
                <div class="d-flex align-items-center">
                  <div style="width: 16px; height: 16px; background-color: ${COMPARISON_CONFIG.colors.unknown}; border-radius: 3px; margin-right: 6px;"></div>
                  <span>Unknown/Potential Sites (${unknownSites.length})</span>
                </div>
              </div>
            </div>
          </div>
          
          <div class="row mb-3">
            <div class="col-lg-6 mb-4">
              <div id="summary-visualization" style="height: 400px;"></div>
            </div>
            <div class="col-lg-6 mb-4">
              <div id="distribution-visualization" style="height: 400px;"></div>
            </div>
          </div>
          
          <div class="row mb-3">
            <div class="col-12">
              <div class="mb-3">
                <label for="metric-selector" class="form-label">Select metric to visualize:</label>
                <select class="form-select" id="metric-selector">
                  ${COMPARISON_CONFIG.metrics.map((metric, index) => 
                    `<option value="${metric.id}" ${index === 0 ? 'selected' : ''}>${metric.name}</option>`
                  ).join('')}
                </select>
              </div>
            </div>
          </div>
          
          <div class="row">
            <div class="col-lg-6 mb-4">
              <div id="boxplot-visualization" style="height: 350px;"></div>
            </div>
            <div class="col-lg-6 mb-4">
              <div id="histogram-visualization" style="height: 350px;"></div>
            </div>
          </div>
          
          <div class="row">
            <div class="col-12 mb-4">
              <div id="siteType-comparison" style="height: 300px;"></div>
            </div>
          </div>
          
          <div class="row">
            <div class="col-12 mb-4">
              <div id="scatter-visualization" style="height: 400px;"></div>
            </div>
          </div>
          
          <div class="row">
            <div class="col-12">
              <p class="text-muted mb-0">
                <em>These visualizations compare the structural and biochemical properties of known phosphorylation sites
                versus potential (not yet verified) sites. Significant differences may indicate structural or
                sequence features that are characteristic of functional phosphorylation sites.</em>
              </p>
            </div>
          </div>
        </div>
      </div>
    `;
    
    // Create the visualizations
    createSummaryVisualization('summary-visualization', knownSites, unknownSites);
    createDistributionVisualization('distribution-visualization', knownSites, unknownSites);
    
    // Create the initial boxplot and histogram with the first metric
    const initialMetric = COMPARISON_CONFIG.metrics[0].id;
    createBoxPlotVisualization('boxplot-visualization', knownSites, unknownSites, initialMetric);
    createHistogramVisualization('histogram-visualization', knownSites, unknownSites, initialMetric);
    
    // Create site type comparison visualization
    createSiteTypeComparison('siteType-comparison', knownSites, unknownSites);
    
    // Create scatter plot visualization (pLDDT vs Surface Accessibility)
    createScatterVisualization('scatter-visualization', knownSites, unknownSites);
    
    // Set up event listener for metric selector
    const metricSelector = document.getElementById('metric-selector');
    if (metricSelector) {
      metricSelector.addEventListener('change', function() {
        const selectedMetric = this.value;
        // Update the visualizations with the new metric
        createBoxPlotVisualization('boxplot-visualization', knownSites, unknownSites, selectedMetric);
        createHistogramVisualization('histogram-visualization', knownSites, unknownSites, selectedMetric);
      });
    }
  }
  
  /**
   * Extract phosphosite data from the table on the page
   * @returns {Array} Array of phosphosite objects
   */
  function extractPhosphositesFromTable() {
    // This function is already implemented in phosphosite-visualization.js
    // but we'll re-implement a simplified version here for our specific needs
    
    console.log('Extracting phosphosites from table for comparison...');
    
    try {
      // Try different table selectors to find phosphosite data
      const tableSelectors = [
        '.phosphosite-table tbody tr',
        'table.table-striped tbody tr',
        '#phosphosite-table tr',
        'table tbody tr'
      ];
      
      let rows = [];
      for (const selector of tableSelectors) {
        const found = document.querySelectorAll(selector);
        if (found && found.length > 0) {
          rows = Array.from(found);
          console.log(`Found ${rows.length} phosphosite rows using selector: ${selector}`);
          break;
        }
      }
      
      if (rows.length === 0) {
        console.warn("No phosphosite table rows found");
        return [];
      }
      
      const phosphosites = [];
      
      // Process each row
      for (const row of rows) {
        try {
          // Initialize phosphosite object with defaults
          const phosphosite = {
            site: '',
            resno: 0,
            siteType: '',
            meanPlddt: 0,
            sitePlddt: 0,
            nearbyCount: 0,
            surfaceAccessibility: 0,
            acidicPercentage: 0,
            basicPercentage: 0,
            isKnown: false
          };
          
          // Get data from data attributes if available
          if (row.hasAttribute('data-site')) {
            phosphosite.site = row.getAttribute('data-site');
            // Extract site type (S, T, Y) from the site attribute
            phosphosite.siteType = phosphosite.site.charAt(0);
          }
          
          if (row.hasAttribute('data-resno')) {
            phosphosite.resno = parseInt(row.getAttribute('data-resno'));
          }
          
          if (row.hasAttribute('data-type')) {
            phosphosite.siteType = row.getAttribute('data-type');
          }
          
          if (row.hasAttribute('data-plddt')) {
            phosphosite.meanPlddt = parseFloat(row.getAttribute('data-plddt'));
          }
          
          if (row.hasAttribute('data-nearby')) {
            phosphosite.nearbyCount = parseInt(row.getAttribute('data-nearby'));
          }
          
          if (row.hasAttribute('data-surface')) {
            phosphosite.surfaceAccessibility = parseFloat(row.getAttribute('data-surface'));
          }
          
          if (row.hasAttribute('data-acidic')) {
            phosphosite.acidicPercentage = parseFloat(row.getAttribute('data-acidic'));
          }
          
          if (row.hasAttribute('data-basic')) {
            phosphosite.basicPercentage = parseFloat(row.getAttribute('data-basic'));
          }
          
          if (row.hasAttribute('data-known')) {
            phosphosite.isKnown = row.getAttribute('data-known') === 'true';
          }
          
          // Extract from table cells as fallback
          const cells = row.querySelectorAll('td');
          if (cells.length >= 3 && !phosphosite.site) {
            // Extract site info from first cell if not already set
            const siteEl = cells[0].querySelector('a') || cells[0];
            const siteText = siteEl.textContent.trim();
            const match = siteText.match(/([STY])(\d+)/);
            
            if (match) {
              phosphosite.site = siteText;
              phosphosite.siteType = match[1];
              phosphosite.resno = parseInt(match[2]);
            }
          }
          
          // Look for special site_plddt field (new field we'll try to extract)
          if (cells.length >= 4 && cells[3]) {
            try {
              const sitePlddt = parseFloat(cells[3].textContent.trim());
              if (!isNaN(sitePlddt)) {
                phosphosite.sitePlddt = sitePlddt;
              }
            } catch (e) {
              // Ignore parsing errors
            }
          }
          
          // If isKnown is not set from data attribute, check "Known" column
          if (phosphosite.isKnown === false) {
            // Look for "Yes" in known column (usually 6th column)
            if (cells.length >= 7 && cells[6]) {
              phosphosite.isKnown = cells[6].textContent.trim() === 'Yes';
            }
            
            // Alternative: check for "Known" in column header and corresponding cell
            for (let i = 0; i < cells.length; i++) {
              if (cells[i] && cells[i].textContent.trim() === 'Yes') {
                // Get the corresponding header cell to see if it's about "Known"
                const headerRow = row.closest('table').querySelector('thead tr');
                if (headerRow) {
                  const headerCells = headerRow.querySelectorAll('th');
                  if (headerCells.length > i) {
                    const headerText = headerCells[i].textContent.toLowerCase();
                    if (headerText.includes('known')) {
                      phosphosite.isKnown = true;
                      break;
                    }
                  }
                }
              }
            }
          }
          
          // Only add valid sites
          if (phosphosite.site && phosphosite.siteType) {
            phosphosites.push(phosphosite);
          }
        } catch (error) {
          console.error("Error processing row:", error);
        }
      }
      
      return phosphosites;
    } catch (error) {
      console.error("Error extracting phosphosites:", error);
      return [];
    }
  }
  
  /**
   * Create a summary visualization with average metrics for known vs unknown sites
   * @param {string} containerId - ID of the container element
   * @param {Array} knownSites - Array of known phosphosite objects
   * @param {Array} unknownSites - Array of unknown phosphosite objects
   */
  function createSummaryVisualization(containerId, knownSites, unknownSites) {
    const container = document.getElementById(containerId);
    if (!container) return;
    
    // Calculate averages for each metric
    const metrics = [
      { id: 'plddt', name: 'Mean pLDDT' },
      { id: 'nearby', name: 'Nearby Residues' },
      { id: 'surface', name: 'Surface Access. (%)' }
    ];
    
    const summaryData = metrics.map(metric => {
      // Get property name based on metric id
      let propName;
      switch (metric.id) {
        case 'plddt': propName = 'meanPlddt'; break;
        case 'nearby': propName = 'nearbyCount'; break;
        case 'surface': propName = 'surfaceAccessibility'; break;
        default: propName = metric.id;
      }
      
      // Calculate averages, removing NaN/null values
      const knownValues = knownSites.map(site => site[propName]).filter(v => !isNaN(v) && v !== null);
      const unknownValues = unknownSites.map(site => site[propName]).filter(v => !isNaN(v) && v !== null);
      
      const knownAvg = knownValues.length > 0 ? knownValues.reduce((a, b) => a + b, 0) / knownValues.length : 0;
      const unknownAvg = unknownValues.length > 0 ? unknownValues.reduce((a, b) => a + b, 0) / unknownValues.length : 0;
      
      return {
        metric: metric.name,
        known: knownAvg,
        unknown: unknownAvg
      };
    });
    
    // Create the summary chart using Chart.js
    const canvas = document.createElement('canvas');
    container.innerHTML = '';
    container.appendChild(canvas);
    
    new Chart(canvas, {
      type: 'bar',
      data: {
        labels: summaryData.map(d => d.metric),
        datasets: [
          {
            label: 'Known Sites',
            data: summaryData.map(d => d.known),
            backgroundColor: COMPARISON_CONFIG.colors.known,
            borderColor: COMPARISON_CONFIG.colors.known,
            borderWidth: 1
          },
          {
            label: 'Unknown Sites',
            data: summaryData.map(d => d.unknown),
            backgroundColor: COMPARISON_CONFIG.colors.unknown,
            borderColor: COMPARISON_CONFIG.colors.unknown,
            borderWidth: 1
          }
        ]
      },
      options: {
        responsive: true,
        maintainAspectRatio: false,
        plugins: {
          title: {
            display: true,
            text: 'Average Metrics Comparison',
            font: {
              size: 16
            }
          },
          legend: {
            position: 'bottom'
          },
          tooltip: {
            callbacks: {
              label: function(context) {
                const value = context.raw.toFixed(1);
                return `${context.dataset.label}: ${value}`;
              }
            }
          }
        },
        scales: {
          y: {
            beginAtZero: true,
            title: {
              display: true,
              text: 'Average Value'
            }
          },
          x: {
            title: {
              display: true,
              text: 'Metric'
            }
          }
        }
      }
    });
  }
  
  /**
   * Create a distribution visualization showing site types (S/T/Y) for known vs unknown sites
   * @param {string} containerId - ID of the container element
   * @param {Array} knownSites - Array of known phosphosite objects
   * @param {Array} unknownSites - Array of unknown phosphosite objects
   */
  function createDistributionVisualization(containerId, knownSites, unknownSites) {
    const container = document.getElementById(containerId);
    if (!container) return;
    
    // Count site types
    const siteTypes = ['S', 'T', 'Y'];
    
    const knownCounts = siteTypes.map(type => 
      knownSites.filter(site => site.siteType === type).length
    );
    
    const unknownCounts = siteTypes.map(type => 
      unknownSites.filter(site => site.siteType === type).length
    );
    
    // Calculate percentages
    const knownTotal = knownCounts.reduce((a, b) => a + b, 0);
    const unknownTotal = unknownCounts.reduce((a, b) => a + b, 0);
    
    const knownPercentages = knownCounts.map(count => (count / knownTotal) * 100);
    const unknownPercentages = unknownCounts.map(count => (count / unknownTotal) * 100);
    
    // Create the pie charts using Chart.js
    container.innerHTML = `
      <div class="row">
        <div class="col-6">
          <canvas id="known-distribution-chart"></canvas>
        </div>
        <div class="col-6">
          <canvas id="unknown-distribution-chart"></canvas>
        </div>
      </div>
    `;
    
    // Colors for site types
    const typeColors = {
      'S': '#4caf50',  // Green for Serine
      'T': '#2196f3',  // Blue for Threonine
      'Y': '#ff9800'   // Orange for Tyrosine
    };
    
    // Create known sites pie chart
    const knownCanvas = document.getElementById('known-distribution-chart');
    new Chart(knownCanvas, {
      type: 'pie',
      data: {
        labels: siteTypes.map(type => `${type} (${knownCounts[siteTypes.indexOf(type)]})`),
        datasets: [{
          data: knownPercentages,
          backgroundColor: siteTypes.map(type => typeColors[type]),
          borderWidth: 1
        }]
      },
      options: {
        responsive: true,
        maintainAspectRatio: false,
        plugins: {
          title: {
            display: true,
            text: 'Known Sites Distribution',
            font: {
              size: 16
            }
          },
          legend: {
            position: 'bottom'
          },
          tooltip: {
            callbacks: {
              label: function(context) {
                const value = context.raw.toFixed(1);
                return `${context.label}: ${value}%`;
              }
            }
          }
        }
      }
    });
    
    // Create unknown sites pie chart
    const unknownCanvas = document.getElementById('unknown-distribution-chart');
    new Chart(unknownCanvas, {
      type: 'pie',
      data: {
        labels: siteTypes.map(type => `${type} (${unknownCounts[siteTypes.indexOf(type)]})`),
        datasets: [{
          data: unknownPercentages,
          backgroundColor: siteTypes.map(type => typeColors[type]),
          borderWidth: 1
        }]
      },
      options: {
        responsive: true,
        maintainAspectRatio: false,
        plugins: {
          title: {
            display: true,
            text: 'Unknown Sites Distribution',
            font: {
              size: 16
            }
          },
          legend: {
            position: 'bottom'
          },
          tooltip: {
            callbacks: {
              label: function(context) {
                const value = context.raw.toFixed(1);
                return `${context.label}: ${value}%`;
              }
            }
          }
        }
      }
    });
  }
  
  /**
   * Create a box plot visualization comparing a metric between known and unknown sites
   * @param {string} containerId - ID of the container element
   * @param {Array} knownSites - Array of known phosphosite objects
   * @param {Array} unknownSites - Array of unknown phosphosite objects
   * @param {string} metricId - ID of the metric to visualize
   */
  function createBoxPlotVisualization(containerId, knownSites, unknownSites, metricId) {
    const container = document.getElementById(containerId);
    if (!container) return;
    
    // Map metric ID to property name
    let propName;
    switch (metricId) {
      case 'plddt': propName = 'meanPlddt'; break;
      case 'site_plddt': propName = 'sitePlddt'; break;
      case 'nearby': propName = 'nearbyCount'; break;
      case 'surface': propName = 'surfaceAccessibility'; break;
      case 'acidic': propName = 'acidicPercentage'; break;
      case 'basic': propName = 'basicPercentage'; break;
      default: propName = metricId;
    }
    
    // Get the metric configuration
    const metricConfig = COMPARISON_CONFIG.metrics.find(m => m.id === metricId) || {
      name: metricId.charAt(0).toUpperCase() + metricId.slice(1),
      format: (value) => value.toFixed(1)
    };
    
    // Extract values, filtering out NaN/null values
    const knownValues = knownSites.map(site => site[propName]).filter(v => !isNaN(v) && v !== null);
    const unknownValues = unknownSites.map(site => site[propName]).filter(v => !isNaN(v) && v !== null);
    
    // If not enough data, show warning
    if (knownValues.length < 5 || unknownValues.length < 5) {
      container.innerHTML = `<div class="alert alert-warning">
        Not enough data to create box plot for ${metricConfig.name}. Need at least 5 values in each group.
      </div>`;
      return;
    }
    
    // Calculate box plot statistics for known sites
    const knownStats = calculateBoxPlotStats(knownValues);
    
    // Calculate box plot statistics for unknown sites
    const unknownStats = calculateBoxPlotStats(unknownValues);
    
    // Create canvas for the chart
    container.innerHTML = '';
    const canvas = document.createElement('canvas');
    container.appendChild(canvas);
    
    // Create the box plot
    new Chart(canvas, {
      type: 'bar',
      data: {
        labels: ['Known Sites', 'Unknown Sites'],
        datasets: [
          // Min to Max lines
          {
            label: 'Min/Max Range',
            data: [
              { x: 'Known Sites', y: knownStats.min, y1: knownStats.max },
              { x: 'Unknown Sites', y: unknownStats.min, y1: unknownStats.max }
            ],
            backgroundColor: 'rgba(0, 0, 0, 0)',
            borderColor: '#666',
            borderWidth: 1,
            barPercentage: 0.15,
            type: 'rangeBar'
          },
          // Q1 to Q3 boxes
          {
            label: 'Interquartile Range',
            data: [
              { x: 'Known Sites', y: knownStats.q1, y1: knownStats.q3 },
              { x: 'Unknown Sites', y: unknownStats.q1, y1: unknownStats.q3 }
            ],
            backgroundColor: [COMPARISON_CONFIG.colors.known, COMPARISON_CONFIG.colors.unknown],
            borderColor: '#333',
            borderWidth: 1,
            barPercentage: 0.4,
            type: 'rangeBar'
          },
          // Median lines
          {
            label: 'Median',
            data: [
              { x: 'Known Sites', y: knownStats.median, y1: knownStats.median },
              { x: 'Unknown Sites', y: unknownStats.median, y1: unknownStats.median }
            ],
            backgroundColor: '#000',
            borderColor: '#000',
            borderWidth: 2,
            barPercentage: 0.5,
            type: 'rangeBar'
          }
        ]
      },
      options: {
        responsive: true,
        maintainAspectRatio: false,
        plugins: {
          title: {
            display: true,
            text: `Box Plot: ${metricConfig.name}`,
            font: {
              size: 16
            }
          },
          legend: {
            display: false
          },
          tooltip: {
            callbacks: {
              label: function(context) {
                const datasetIndex = context.datasetIndex;
                const dataIndex = context.dataIndex;
                
                if (datasetIndex === 0) { // Min/Max
                  const min = context.dataset.data[dataIndex].y;
                  const max = context.dataset.data[dataIndex].y1;
                  return [
                    `Min: ${metricConfig.format(min)}`,
                    `Max: ${metricConfig.format(max)}`
                  ];
                } else if (datasetIndex === 1) { // IQR
                  const q1 = context.dataset.data[dataIndex].y;
                  const q3 = context.dataset.data[dataIndex].y1;
                  return [
                    `Q1 (25%): ${metricConfig.format(q1)}`,
                    `Q3 (75%): ${metricConfig.format(q3)}`
                  ];
                } else { // Median
                  const median = context.dataset.data[dataIndex].y;
                  return `Median: ${metricConfig.format(median)}`;
                }
              }
            }
          }
        },
        scales: {
          y: {
            beginAtZero: metricId !== 'plddt' && metricId !== 'site_plddt',
            title: {
              display: true,
              text: metricConfig.name
            }
          }
        }
      }
    });
    
    // Add description below the chart
    const descriptionDiv = document.createElement('div');
    descriptionDiv.className = 'text-muted mt-2 small';
    descriptionDiv.innerHTML = `
      <p><strong>Box Plot:</strong> ${metricConfig.description}</p>
      <p><strong>Known Sites:</strong> Median = ${metricConfig.format(knownStats.median)}, 
         IQR = ${metricConfig.format(knownStats.q1)} - ${metricConfig.format(knownStats.q3)}</p>
      <p><strong>Unknown Sites:</strong> Median = ${metricConfig.format(unknownStats.median)}, 
         IQR = ${metricConfig.format(unknownStats.q1)} - ${metricConfig.format(unknownStats.q3)}</p>
    `;
    container.appendChild(descriptionDiv);
  }
  
  /**
   * Check if Chart.js is loaded and load it if necessary
   * @param {Function} callback - Function to call when Chart.js is loaded
   */
  // Check if Chart.js is loaded and load it if necessary
function ensureChartJsLoaded(callback) {
    if (typeof Chart !== 'undefined') {
      // Chart.js is already loaded
      if (callback) callback();
      return;
    }
    
    console.log('Loading Chart.js dynamically');
    
    // Load Chart.js
    const script = document.createElement('script');
    script.src = 'https://cdn.jsdelivr.net/npm/chart.js@3.9.1/dist/chart.min.js';
    script.onload = function() {
      console.log('Chart.js loaded successfully');
      
      // Load Chart.js range/box-plot plugin
      const rangePlugin = document.createElement('script');
      rangePlugin.src = 'https://cdn.jsdelivr.net/npm/chartjs-chart-box-and-violin-plot@4.0.0/build/index.umd.min.js';
      rangePlugin.onload = function() {
        // After both are loaded, call the callback
        if (callback) callback();
      };
      document.head.appendChild(rangePlugin);
    };
    
    script.onerror = function() {
      console.error('Failed to load Chart.js dynamically');
      alert('Failed to load Chart.js library. Visualizations may not work properly.');
    };
    
    document.head.appendChild(script);
  }
  
  // -- AUTO-INITIALIZE ON PAGE LOAD --
document.addEventListener('DOMContentLoaded', function() {
    // Look for comparison container
    const comparisonContainer = document.getElementById('phosphosite-comparison-container');
    if (comparisonContainer) {
      // Show the container (which may be hidden by default)
      comparisonContainer.style.display = 'block';
      
      // Ensure Chart.js is loaded before initializing
      ensureChartJsLoaded(function() {
        // Initialize comparison visualization
        initPhosphositeComparison('phosphosite-comparison-container');
      });
    } else {
      console.warn('Phosphosite comparison container not found on the page');
    }
    
    // Add toggle buttons if the visualization container exists
    const vizContainer = document.getElementById('phosphosite-visualization-container');
    if (vizContainer) {
      // Create a button group for visualization options
      const buttonGroup = document.createElement('div');
      buttonGroup.className = 'btn-group mb-3';
      buttonGroup.role = 'group';
      buttonGroup.setAttribute('aria-label', 'Visualization Options');
      
      // Create standard visualization button
      const standardButton = document.createElement('button');
      standardButton.type = 'button';
      standardButton.className = 'btn btn-primary active';
      standardButton.textContent = 'Standard View';
      standardButton.onclick = function() {
        // Show standard container
        vizContainer.style.display = 'block';
        
        // Hide comparison container if it exists
        const comparisonContainer = document.getElementById('phosphosite-comparison-container');
        if (comparisonContainer) {
          comparisonContainer.style.display = 'none';
        }
        
        // Update active button
        standardButton.classList.add('active');
        comparisonButton.classList.remove('active');
        
        // Trigger window resize to refresh charts
        window.dispatchEvent(new Event('resize'));
      };
      
      // Create comparison visualization button
      const comparisonButton = document.createElement('button');
      comparisonButton.type = 'button';
      comparisonButton.className = 'btn btn-outline-primary';
      comparisonButton.textContent = 'Known vs Unknown Comparison';
      comparisonButton.onclick = function() {
        // Hide standard container
        vizContainer.style.display = 'none';
        
        // Show comparison container
        const comparisonContainer = document.getElementById('phosphosite-comparison-container');
        if (comparisonContainer) {
          comparisonContainer.style.display = 'block';
          
          // Initialize comparison visualizations if needed
          if (comparisonContainer.innerHTML === '' || comparisonContainer.innerHTML.includes('Comparison Visualization Container')) {
            ensureChartJsLoaded(function() {
              initPhosphositeComparison('phosphosite-comparison-container');
            });
          }
        } else {
          console.warn('Comparison container not found');
        }
        
        // Update active button
        standardButton.classList.remove('active');
        comparisonButton.classList.add('active');
        
        // Trigger window resize to refresh charts
        window.dispatchEvent(new Event('resize'));
      };
      
      // Add buttons to the group
      buttonGroup.appendChild(standardButton);
      buttonGroup.appendChild(comparisonButton);
      
      // Add the button group to the page
      const phosphositeTable = document.querySelector('.phosphosite-table');
      if (phosphositeTable) {
        // If table exists, add buttons before it
        const tableCard = phosphositeTable.closest('.card');
        if (tableCard) {
          tableCard.parentNode.insertBefore(buttonGroup, tableCard);
        } else {
          // Fallback - add buttons before the visualization container
          vizContainer.parentNode.insertBefore(buttonGroup, vizContainer);
        }
      } else {
        // Fallback - add buttons before the visualization container
        vizContainer.parentNode.insertBefore(buttonGroup, vizContainer);
      }
    }
  });
  
  // Export the initialization function for external use
  window.initPhosphositeComparison = initPhosphositeComparison;
  
  /**
   * Calculate statistics for a box plot
   * @param {Array} values - Array of numeric values
   * @returns {Object} Object with min, q1, median, q3, max
   */
  function calculateBoxPlotStats(values) {
    const sortedValues = [...values].sort((a, b) => a - b);
    const n = sortedValues.length;
    
    // Get quartiles
    const q1Index = Math.floor(n * 0.25);
    const medianIndex = Math.floor(n * 0.5);
    const q3Index = Math.floor(n * 0.75);
    
    return {
      min: sortedValues[0],
      q1: sortedValues[q1Index],
      median: sortedValues[medianIndex],
      q3: sortedValues[q3Index],
      max: sortedValues[n - 1]
    };
  }
  
  /**
   * Create a histogram visualization comparing a metric between known and unknown sites
   * @param {string} containerId - ID of the container element
   * @param {Array} knownSites - Array of known phosphosite objects
   * @param {Array} unknownSites - Array of unknown phosphosite objects
   * @param {string} metricId - ID of the metric to visualize
   */
  function createHistogramVisualization(containerId, knownSites, unknownSites, metricId) {
    const container = document.getElementById(containerId);
    if (!container) return;
    
    // Map metric ID to property name
    let propName;
    switch (metricId) {
      case 'plddt': propName = 'meanPlddt'; break;
      case 'site_plddt': propName = 'sitePlddt'; break;
      case 'nearby': propName = 'nearbyCount'; break;
      case 'surface': propName = 'surfaceAccessibility'; break;
      case 'acidic': propName = 'acidicPercentage'; break;
      case 'basic': propName = 'basicPercentage'; break;
      default: propName = metricId;
    }
    
    // Get the metric configuration
    const metricConfig = COMPARISON_CONFIG.metrics.find(m => m.id === metricId) || {
      name: metricId.charAt(0).toUpperCase() + metricId.slice(1),
      format: (value) => value.toFixed(1)
    };
    
    // Extract values, filtering out NaN/null values
    const knownValues = knownSites.map(site => site[propName]).filter(v => !isNaN(v) && v !== null);
    const unknownValues = unknownSites.map(site => site[propName]).filter(v => !isNaN(v) && v !== null);
    
    // If not enough data, show warning
    if (knownValues.length === 0 || unknownValues.length === 0) {
      container.innerHTML = `<div class="alert alert-warning">
        Not enough data to create histogram for ${metricConfig.name}.
      </div>`;
      return;
    }
    
    // Determine bin range and count
    const allValues = [...knownValues, ...unknownValues];
    const min = Math.min(...allValues);
    const max = Math.max(...allValues);
    const range = max - min;
    
    // Determine number of bins - between 6 and 12 depending on data size
    const minBins = 6;
    const maxBins = 12;
    const dataSize = allValues.length;
    const binCount = Math.min(maxBins, Math.max(minBins, Math.floor(Math.sqrt(dataSize))));
    
    // Calculate bin width
    const binWidth = range / binCount;
    
    // Create histogram bins
    const bins = Array(binCount).fill(0).map((_, i) => ({
      start: min + i * binWidth,
      end: min + (i + 1) * binWidth,
      knownCount: 0,
      unknownCount: 0
    }));
    
    // Fill the bins
    knownValues.forEach(value => {
      const binIndex = Math.min(binCount - 1, Math.floor((value - min) / binWidth));
      bins[binIndex].knownCount++;
    });
    
    unknownValues.forEach(value => {
      const binIndex = Math.min(binCount - 1, Math.floor((value - min) / binWidth));
      bins[binIndex].unknownCount++;
    });
    
    // Convert to percentages for better comparison
    bins.forEach(bin => {
      bin.knownPercent = (bin.knownCount / knownValues.length) * 100;
      bin.unknownPercent = (bin.unknownCount / unknownValues.length) * 100;
    });
    
    // Create canvas for the chart
    container.innerHTML = '';
    const canvas = document.createElement('canvas');
    container.appendChild(canvas);
    
    // Create the histogram
    new Chart(canvas, {
      type: 'bar',
      data: {
        labels: bins.map(bin => `${metricConfig.format(bin.start)} - ${metricConfig.format(bin.end)}`),
        datasets: [
          {
            label: 'Known Sites (%)',
            data: bins.map(bin => bin.knownPercent),
            backgroundColor: COMPARISON_CONFIG.colors.known,
            borderColor: COMPARISON_CONFIG.colors.known,
            borderWidth: 1
          },
          {
            label: 'Unknown Sites (%)',
            data: bins.map(bin => bin.unknownPercent),
            backgroundColor: COMPARISON_CONFIG.colors.unknown,
            borderColor: COMPARISON_CONFIG.colors.unknown,
            borderWidth: 1
          }
        ]
      },
      options: {
        responsive: true,
        maintainAspectRatio: false,
        plugins: {
          title: {
            display: true,
            text: `Distribution: ${metricConfig.name}`,
            font: {
              size: 16
            }
          },
          legend: {
            position: 'bottom'
          },
          tooltip: {
            callbacks: {
              label: function(context) {
                const datasetIndex = context.datasetIndex;
                const binIndex = context.dataIndex;
                const bin = bins[binIndex];
                
                if (datasetIndex === 0) { // Known sites
                  return [
                    `Known Sites: ${bin.knownPercent.toFixed(1)}%`,
                    `Count: ${bin.knownCount} of ${knownValues.length}`
                  ];
                } else { // Unknown sites
                  return [
                    `Unknown Sites: ${bin.unknownPercent.toFixed(1)}%`,
                    `Count: ${bin.unknownCount} of ${unknownValues.length}`
                  ];
                }
              }
            }
          }
        },
        scales: {
          y: {
            beginAtZero: true,
            title: {
              display: true,
              text: 'Percentage (%)'
            }
          },
          x: {
            title: {
              display: true,
              text: metricConfig.name
            }
          }
        }
      }
    });
  }
  
  /**
   * Create a visualization comparing site types between known and unknown sites
   * @param {string} containerId - ID of the container element
   * @param {Array} knownSites - Array of known phosphosite objects
   * @param {Array} unknownSites - Array of unknown phosphosite objects
   */
  function createSiteTypeComparison(containerId, knownSites, unknownSites) {
    const container = document.getElementById(containerId);
    if (!container) return;
    
    // Count site types by Serine, Threonine, Tyrosine
    const siteTypes = ['S', 'T', 'Y'];
    const typeLabels = ['Serine (S)', 'Threonine (T)', 'Tyrosine (Y)'];
    
    // Count for known sites
    const knownTypeCounts = siteTypes.map(type => 
      knownSites.filter(site => site.siteType === type).length
    );
    
    // Count for unknown sites
    const unknownTypeCounts = siteTypes.map(type => 
      unknownSites.filter(site => site.siteType === type).length
    );
    
    // Calculate percentages
    const knownTotal = knownSites.length;
    const unknownTotal = unknownSites.length;
    
    const knownPercents = knownTypeCounts.map(count => (count / knownTotal) * 100);
    const unknownPercents = unknownTypeCounts.map(count => (count / unknownTotal) * 100);
    
    // Create canvas for the chart
    container.innerHTML = '';
    const canvas = document.createElement('canvas');
    container.appendChild(canvas);
    
    // Create the comparison chart
    new Chart(canvas, {
      type: 'bar',
      data: {
        labels: typeLabels,
        datasets: [
          {
            label: 'Known Sites (%)',
            data: knownPercents,
            backgroundColor: COMPARISON_CONFIG.colors.known,
            borderColor: COMPARISON_CONFIG.colors.known,
            borderWidth: 1
          },
          {
            label: 'Unknown Sites (%)',
            data: unknownPercents,
            backgroundColor: COMPARISON_CONFIG.colors.unknown,
            borderColor: COMPARISON_CONFIG.colors.unknown,
            borderWidth: 1
          }
        ]
      },
      options: {
        responsive: true,
        maintainAspectRatio: false,
        plugins: {
          title: {
            display: true,
            text: 'Site Type Distribution Comparison',
            font: {
              size: 16
            }
          },
          legend: {
            position: 'bottom'
          },
          tooltip: {
            callbacks: {
              label: function(context) {
                const datasetIndex = context.datasetIndex;
                const index = context.dataIndex;
                const percent = context.raw.toFixed(1);
                
                if (datasetIndex === 0) { // Known sites
                  const count = knownTypeCounts[index];
                  return [
                    `Known Sites: ${percent}%`,
                    `Count: ${count} of ${knownTotal}`
                  ];
                } else { // Unknown sites
                  const count = unknownTypeCounts[index];
                  return [
                    `Unknown Sites: ${percent}%`,
                    `Count: ${count} of ${unknownTotal}`
                  ];
                }
              }
            }
          }
        },
        scales: {
          y: {
            beginAtZero: true,
            title: {
              display: true,
              text: 'Percentage (%)'
            }
          }
        }
      }
    });
  }
  
  /**
   * Create a scatter plot visualization comparing two metrics
   * @param {string} containerId - ID of the container element
   * @param {Array} knownSites - Array of known phosphosite objects
   * @param {Array} unknownSites - Array of unknown phosphosite objects
   */
  function createScatterVisualization(containerId, knownSites, unknownSites) {
    const container = document.getElementById(containerId);
    if (!container) return;
    
    // We'll create a scatter plot of pLDDT vs Surface Accessibility
    const xMetric = 'meanPlddt';
    const yMetric = 'surfaceAccessibility';
    
    // Extract values for known sites
    const knownData = knownSites.map(site => ({
      x: site[xMetric] || 0,
      y: site[yMetric] || 0,
      site: site.site,
      siteType: site.siteType
    })).filter(point => !isNaN(point.x) && !isNaN(point.y));
    
    // Extract values for unknown sites
    const unknownData = unknownSites.map(site => ({
      x: site[xMetric] || 0,
      y: site[yMetric] || 0,
      site: site.site,
      siteType: site.siteType
    })).filter(point => !isNaN(point.x) && !isNaN(point.y));
    
    // If not enough data, show warning
    if (knownData.length < 5 || unknownData.length < 5) {
      container.innerHTML = `<div class="alert alert-warning">
        Not enough data to create scatter plot. Need at least 5 values in each group.
      </div>`;
      return;
    }
    
    // Color mapping for site types
    const siteTypeColors = {
      'S': 'rgba(76, 175, 80, 0.7)',  // Green for Serine
      'T': 'rgba(33, 150, 243, 0.7)',  // Blue for Threonine
      'Y': 'rgba(255, 152, 0, 0.7)'    // Orange for Tyrosine
    };
    
    // Create canvas for the chart
    container.innerHTML = '';
    const canvas = document.createElement('canvas');
    container.appendChild(canvas);
    
    // Create the scatter plot
    new Chart(canvas, {
      type: 'scatter',
      data: {
        datasets: [
          {
            label: 'Known Sites',
            data: knownData,
            backgroundColor: knownData.map(point => siteTypeColors[point.siteType] || 'rgba(0, 0, 0, 0.7)'),
            borderColor: 'rgba(0, 0, 0, 0.5)',
            borderWidth: 1,
            pointRadius: 5,
            pointHoverRadius: 7
          },
          {
            label: 'Unknown Sites',
            data: unknownData,
            backgroundColor: function(context) {
              const index = context.dataIndex;
              const siteType = unknownData[index].siteType;
              // Make unknown sites partially transparent
              const baseColor = siteTypeColors[siteType] || 'rgba(0, 0, 0, 0.7)';
              return baseColor.replace('0.7', '0.3');
            },
            borderColor: 'rgba(0, 0, 0, 0.3)',
            borderWidth: 1,
            pointRadius: 4,
            pointHoverRadius: 6
          }
        ]
      },
      options: {
        responsive: true,
        maintainAspectRatio: false,
        plugins: {
          title: {
            display: true,
            text: 'pLDDT Score vs Surface Accessibility',
            font: {
              size: 16
            }
          },
          legend: {
            position: 'bottom'
          },
          tooltip: {
            callbacks: {
              label: function(context) {
                const datasetIndex = context.datasetIndex;
                const index = context.dataIndex;
                
                if (datasetIndex === 0) { // Known sites
                  const point = knownData[index];
                  return [
                    `Site: ${point.site}`,
                    `pLDDT: ${point.x.toFixed(1)}`,
                    `Surface Acc: ${point.y.toFixed(1)}%`
                  ];
                } else { // Unknown sites
                  const point = unknownData[index];
                  return [
                    `Site: ${point.site}`,
                    `pLDDT: ${point.x.toFixed(1)}`,
                    `Surface Acc: ${point.y.toFixed(1)}%`
                  ];
                }
              }
            }
          }
        },
        scales: {
          y: {
            title: {
              display: true,
              text: 'Surface Accessibility (%)'
            }
          },
          x: {
            title: {
              display: true,
              text: 'pLDDT Score'
            }
          }
        }
      }
    });
    
    // Add legend for site types
    const legendDiv = document.createElement('div');
    legendDiv.className = 'd-flex justify-content-center mt-3';
    legendDiv.innerHTML = `
      <div class="d-flex align-items-center me-4">
        <div style="width: 12px; height: 12px; background-color: ${siteTypeColors['S']}; border-radius: 50%; margin-right: 5px;"></div>
        <span>Serine (S)</span>
      </div>
      <div class="d-flex align-items-center me-4">
        <div style="width: 12px; height: 12px; background-color: ${siteTypeColors['T']}; border-radius: 50%; margin-right: 5px;"></div>
        <span>Threonine (T)</span>
      </div>
      <div class="d-flex align-items-center">
        <div style="width: 12px; height: 12px; background-color: ${siteTypeColors['Y']}; border-radius: 50%; margin-right: 5px;"></div>
        <span>Tyrosine (Y)</span>
      </div>
    `;
    container.appendChild(legendDiv);
    
    // Add description
    const descriptionDiv = document.createElement('div');
    descriptionDiv.className = 'text-muted mt-2 small';
    descriptionDiv.innerHTML = `
      <p>This scatter plot shows the relationship between structure quality (pLDDT score) 
      and surface accessibility for known vs unknown phosphorylation sites. Known functional sites 
      (solid colors) often show different structural patterns compared to unverified sites (transparent).</p>
    `;
    container.appendChild(descriptionDiv);
  }