// Phosphosite UI Enhancements - Show both visualizations
// This script modifies the display logic to show both standard and comparison visualizations

document.addEventListener('DOMContentLoaded', function() {
    // Get visualization containers
    const standardContainer = document.getElementById('phosphosite-visualization-container');
    const comparisonContainer = document.getElementById('phosphosite-comparison-container');
    
    // Function to initialize both visualizations when data is available
    function initializeVisualizations() {
      // First check if we have valid content to show
      if (comparisonContainer && standardContainer) {
        // Set a better title for the comparison section
        if (!document.getElementById('comparison-section-title')) {
          const comparisonTitle = document.createElement('h3');
          comparisonTitle.id = 'comparison-section-title';
          comparisonTitle.className = 'mt-5 mb-3';
          comparisonTitle.innerHTML = 'Known vs Unknown Phosphosite Comparison';
          
          // Insert before the comparison container
          if (comparisonContainer.parentNode) {
            comparisonContainer.parentNode.insertBefore(comparisonTitle, comparisonContainer);
          }
        }
        
        // Show both containers
        standardContainer.style.display = 'block';
        comparisonContainer.style.display = 'block';
      }
    }
    
    // Initialize both visualizations on page load
    initializeVisualizations();
    
    // Remove the toggle buttons since we'll show both views
    const toggleButtons = document.querySelectorAll('.view-selector, .btn-group button[data-view]');
    toggleButtons.forEach(button => {
      button.style.display = 'none';
    });
    
    // Check if we need to initialize the visualizations
    const hasPfx = document.querySelector('.phosphosite-table') || 
                  document.getElementById('phosphosite-table') ||
                  document.querySelector('.results-container');
                  
    if (hasPfx) {
      console.log('Found phosphosite data, initializing visualizations');
      
      // Make sure the chart.js library is loaded
      if (typeof Chart === 'undefined') {
        console.log('Chart.js not loaded, loading it now');
        loadChartJs(() => {
          // Once loaded, initialize both visualizations
          if (typeof window.createPhosphositeVisualization === 'function') {
            window.createPhosphositeVisualization('phosphosite-visualization-container');
          }
          
          if (typeof window.initPhosphositeComparison === 'function') {
            window.initPhosphositeComparison('phosphosite-comparison-container');
          } else {
            // Try to load the comparison script
            loadComparisonScript();
          }
        });
      } else {
        // Chart.js is already loaded
        console.log('Chart.js already loaded');
        
        // Initialize standard visualization if not already
        if (standardContainer && standardContainer.innerHTML.trim() === '') {
          if (typeof window.createPhosphositeVisualization === 'function') {
            window.createPhosphositeVisualization('phosphosite-visualization-container');
          }
        }
        
        // Initialize comparison visualization if not already
        if (comparisonContainer && 
            (comparisonContainer.innerHTML.trim() === '' || 
             comparisonContainer.innerHTML.includes('Loading comparison'))) {
          if (typeof window.initPhosphositeComparison === 'function') {
            window.initPhosphositeComparison('phosphosite-comparison-container');
          } else {
            // Try to load the comparison script
            loadComparisonScript();
          }
        }
      }
    }
    
    function loadChartJs(callback) {
      const script = document.createElement('script');
      script.src = 'https://cdn.jsdelivr.net/npm/chart.js@3.9.1/dist/chart.min.js';
      script.onload = function() {
        console.log('Chart.js loaded successfully');
        
        // Load the annotation plugin
        const annotationPlugin = document.createElement('script');
        annotationPlugin.src = 'https://cdn.jsdelivr.net/npm/chartjs-plugin-annotation@2.0.1/dist/chartjs-plugin-annotation.min.js';
        annotationPlugin.onload = function() {
          if (callback) callback();
        };
        document.head.appendChild(annotationPlugin);
      };
      
      document.head.appendChild(script);
    }
    
    function loadComparisonScript() {
      const script = document.createElement('script');
      script.src = "/static/js/phosphosite-comparison.js";
      script.onload = function() {
        console.log('Loaded phosphosite-comparison.js');
        if (typeof window.initPhosphositeComparison === 'function') {
          window.initPhosphositeComparison('phosphosite-comparison-container');
        }
      };
      document.head.appendChild(script);
    }
    
    // Observer to watch for content changes in the containers
    const observer = new MutationObserver(function(mutations) {
      mutations.forEach(function(mutation) {
        // If containers are populated with meaningful content, ensure they're both visible
        if ((standardContainer && standardContainer.innerHTML.length > 100 && 
             !standardContainer.innerHTML.includes('No phosphosite data available')) ||
            (comparisonContainer && comparisonContainer.innerHTML.length > 100 && 
             !comparisonContainer.innerHTML.includes('No phosphosite data available'))) {
          
          // Show both containers
          if (standardContainer) standardContainer.style.display = 'block';
          if (comparisonContainer) comparisonContainer.style.display = 'block';
          
          // Add comparison title if not already present
          if (!document.getElementById('comparison-section-title') && comparisonContainer) {
            const comparisonTitle = document.createElement('h3');
            comparisonTitle.id = 'comparison-section-title';
            comparisonTitle.className = 'mt-5 mb-3';
            comparisonTitle.innerHTML = 'Known vs Unknown Phosphosite Comparison';
            
            // Insert before the comparison container
            if (comparisonContainer.parentNode) {
              comparisonContainer.parentNode.insertBefore(comparisonTitle, comparisonContainer);
            }
          }
        }
      });
    });
    
    // Observe both containers for changes
    if (standardContainer) {
      observer.observe(standardContainer, { childList: true, subtree: true });
    }
    
    if (comparisonContainer) {
      observer.observe(comparisonContainer, { childList: true, subtree: true });
    }
  });