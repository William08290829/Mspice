function toggleMenu() {
  const menu = document.querySelector(".menu-links");
  const icon = document.querySelector(".hamburger-icon");
  menu.classList.toggle("open");
  icon.classList.toggle("open");
}


// SLIDER
const sections = document.querySelectorAll('section[id]');

const scrollActive = () => {
  const scrollDown = window.scrollY

  sections.forEach(current => {
    const sectionHeight = current.offsetHeight;
    const sectionTop = current.offsetTop - window.innerHeight * 0.4;
    const sectionId = current.getAttribute('id');
    let sectionsClass = document.querySelector('.page-slider a[href*=' + sectionId + ']')

    if (sectionsClass === null) {
      sectionsClass = document.querySelector('.page-slider a[id*=' + sectionId + ']')
    }

    if (sectionId === 'desktop-nav') {
      // console.log(scrollDown)
      if (scrollDown <= sectionTop + sectionHeight) {
        sectionsClass.classList.add('paw-active')
      } else {
        sectionsClass.classList.remove('paw-active')
      }
    } else {
      if (scrollDown > sectionTop && scrollDown <= sectionTop + sectionHeight) {
        sectionsClass.classList.add('paw-active')
      } else {
        sectionsClass.classList.remove('paw-active')
      }
    }
  })
}

window.addEventListener('scroll', scrollActive);



// SCROLL ANIMATIONS
const observer = new IntersectionObserver((entries) => {
  entries.forEach((entry) => {
    if (entry.isIntersecting) {
      entry.target.classList.add('show');

      // Trigger typewriter if it's the code block
      if (entry.target.querySelector('.spFile-text') && !typeWriterStarted) {
        typeWriter();
        typeWriterStarted = true;
      }
    }
  });
});

const hiddenElements = document.querySelectorAll('.home__text, .home__pic-container, .intro__container, .about__text, .bread-board-edge, .about__member, .get-started__container, .contact-info-upper-container, section#docs');
hiddenElements.forEach((el) => {
  el.classList.add('hidden');
  observer.observe(el);
});


// TYPEWRITER EFFECT
const codeText = `* Mu circuit

* Netlist
L_head left-ear right-ear 100u
C_left-whisker left-ear node_1 100u
C_left-whisker right-ear node_2 100u
R_mouth left-cheek right-cheek 1k

* Analysis
V1 node_1 0 1
V2 node_2 0 0

.end
`;

let i = 0;
let typeWriterStarted = false;
const speed = 20; /* The speed/duration of the effect in milliseconds */
const resetDelay = 3000; /* Delay before restarting */

function typeWriter() {
  if (i < codeText.length) {
    let char = codeText.charAt(i);
    if (char === '\n') {
      char = '<br>';
    }
    document.querySelector(".spFile-text").innerHTML += char;
    i++;
    setTimeout(typeWriter, speed);
  } else {
    // Reset and loop after delay
    setTimeout(() => {
      document.querySelector(".spFile-text").innerHTML = '';
      i = 0;
      typeWriter();
    }, resetDelay);
  }
}


// ==========================================
// APP LOGIC (app.html)
// ==========================================

const loadingOverlay = document.getElementById('loading-overlay');
const appInterface = document.getElementById('app-interface');
const netlistInput = document.getElementById('netlist-input');
const fileInput = document.getElementById('file-input');
const outputArea = document.getElementById('output-area');

// 1. Preloader Logic
if (loadingOverlay) {
  window.addEventListener('load', () => {
    // Fake delay for effect (1.5s)
    setTimeout(() => {
      loadingOverlay.classList.add('fade-out');
      if (appInterface) appInterface.classList.remove('hidden-app');
    }, 1500);
  });
}

// 2. File Upload Logic
if (fileInput && netlistInput) {
  fileInput.addEventListener('change', (e) => {
    const file = e.target.files[0];
    if (!file) return;

    const reader = new FileReader();
    reader.onload = (event) => {
      netlistInput.value = event.target.result;
    };
    reader.readAsText(file);
  });
}

// 3. Run Simulation Logic
// Global variables to store results for download
let globalPlotImage = null;
let globalCSVData = null;

// 3. Run Simulation Logic
// Global variable for Plotly
let plotlyGraphDiv = null;

async function runSimulation() {
  if (!netlistInput || !outputArea) return;

  const netlist = netlistInput.value;
  if (!netlist.trim()) {
    alert("Please enter or upload a netlist first!");
    return;
  }

  // Reset UI
  document.getElementById('download-plot-btn').classList.add('app-hidden');
  document.getElementById('download-csv-btn').classList.add('app-hidden');
  document.getElementById('data-table').classList.add('app-hidden');
  document.getElementById('signal-controls').classList.add('app-hidden');
  document.getElementById('signal-checkboxes').innerHTML = '';

  const dataPlaceholder = document.getElementById('data-placeholder');
  if (dataPlaceholder) {
    dataPlaceholder.innerHTML = '<div class="spinner" style="width:30px;height:30px;"></div><p>Simulating...</p>';
    dataPlaceholder.classList.remove('app-hidden');
  }

  globalPlotImage = null;
  globalCSVData = null;

  // Show loading state in output
  outputArea.innerHTML = '<div class="result-placeholder"><div class="spinner" style="width:30px;height:30px;"></div><p>Simulating...</p></div>';

  // Button Loading State
  const runBtn = document.querySelector('.btn-app.primary');
  const originalBtnContent = runBtn.innerHTML;

  if (runBtn) {
    runBtn.innerHTML = '<div class="mini-spinner"></div> Running...';
    runBtn.classList.add('btn-loading');
    runBtn.disabled = true;
  }

  try {
    const response = await fetch('/api/simulate', {
      method: 'POST',
      headers: {
        'Content-Type': 'application/json'
      },
      body: JSON.stringify({ netlist: netlist })
    });

    const data = await response.json();

    if (response.ok && data.data) {
      // 1. Render Plotly & Controls
      globalCSVData = data.data;
      renderPlotly(data.data);

      document.getElementById('download-plot-btn').classList.remove('app-hidden');
      document.getElementById('signal-controls').classList.remove('app-hidden');

      // 2. Render Data Table & Setup CSV Download
      renderTable(data.data);

      if (document.getElementById('data-placeholder')) document.getElementById('data-placeholder').classList.add('app-hidden');
      document.getElementById('data-table').classList.remove('app-hidden');
      document.getElementById('download-csv-btn').classList.remove('app-hidden');

    } else if (data.error) {
      outputArea.innerHTML = `<p style="color:red;">Error: ${data.error}</p>`;
    } else {
      outputArea.innerHTML = `<p>${data.message || 'Unknown result'}</p>`;
    }

  } catch (err) {
    console.error(err);
    outputArea.innerHTML = `<p style="color:red;">Network Error: ${err.message}</p>`;
  } finally {
    // Restore Placeholder Text
    const dataPlaceholder = document.getElementById('data-placeholder');
    if (dataPlaceholder) {
      dataPlaceholder.innerHTML = '<p class="placeholder-text">Run a simulation to see raw data here.</p>';
    }

    // Reset Button State
    if (runBtn) {
      runBtn.innerHTML = originalBtnContent;
      runBtn.classList.remove('btn-loading');
      runBtn.disabled = false;

      // Suppress animation until mouse leaves
      runBtn.classList.add('suppress-peek');
      runBtn.addEventListener('mouseleave', () => {
        runBtn.classList.remove('suppress-peek');
      }, { once: true });
    }
  }
}

function renderPlotly(csvData) {
  const outputArea = document.getElementById('output-area');
  outputArea.innerHTML = '';

  const plotDiv = document.createElement('div');
  plotDiv.id = 'plotly-chart';
  plotDiv.style.width = '100%';
  plotDiv.style.height = '100%';
  outputArea.appendChild(plotDiv);
  plotlyGraphDiv = plotDiv;

  const headers = csvData.headers;
  const rows = csvData.rows;

  const columns = headers.map(() => []);
  rows.forEach(row => {
    row.forEach((val, idx) => columns[idx].push(val));
  });

  const xData = columns[0];
  const xLabel = headers[0];

  const traces = [];
  const signals = [];

  for (let i = 1; i < headers.length; i++) {
    const trace = {
      x: xData,
      y: columns[i],
      mode: 'lines',
      name: headers[i],
      visible: true
    };
    traces.push(trace);
    signals.push({ index: i - 1, name: headers[i] });
  }

  const layout = {
    margin: { t: 20, r: 20, b: 50, l: 60 },
    xaxis: { title: xLabel },
    yaxis: { title: 'Value' },
    showlegend: false,
    autosize: true
  };

  const config = { responsive: true, displayModeBar: false };

  Plotly.newPlot(plotDiv, traces, layout, config);

  // --- POPULATE DROPDOWN ---
  const checkboxContainer = document.getElementById('signal-checkboxes');
  checkboxContainer.innerHTML = '';

  const signalBoxes = [];

  // 1. Select All Option
  const selectAllLabel = document.createElement('label');
  selectAllLabel.style.display = 'flex';
  selectAllLabel.style.alignItems = 'center';
  selectAllLabel.style.cursor = 'pointer';
  selectAllLabel.style.padding = '4px 0';
  selectAllLabel.style.borderBottom = '1px solid #eee';
  selectAllLabel.style.marginBottom = '5px';
  selectAllLabel.style.fontWeight = 'bold';
  selectAllLabel.style.fontSize = '0.9em';

  const selectAllBox = document.createElement('input');
  selectAllBox.type = 'checkbox';
  selectAllBox.checked = true;
  selectAllBox.style.marginRight = '8px';

  selectAllLabel.appendChild(selectAllBox);
  selectAllLabel.appendChild(document.createTextNode("Select All"));
  checkboxContainer.appendChild(selectAllLabel);

  // 2. Individual Options
  signals.forEach((sig, idx) => {
    const label = document.createElement('label');
    label.style.display = 'flex';
    label.style.alignItems = 'center';
    label.style.cursor = 'pointer';
    label.style.padding = '2px 0';
    label.style.fontSize = '0.9em';

    const box = document.createElement('input');
    box.type = 'checkbox';
    box.checked = true;
    box.style.marginRight = '8px';

    // Individual Change
    box.addEventListener('change', (e) => {
      const isChecked = e.target.checked;
      Plotly.restyle(plotDiv, { visible: isChecked }, [sig.index]);

      // Update Select All State
      if (!isChecked) {
        selectAllBox.checked = false;
      } else {
        // Check if all are now checked
        const allChecked = signalBoxes.every(b => b.checked);
        if (allChecked) selectAllBox.checked = true;
      }
    });

    signalBoxes.push(box);
    label.appendChild(box);
    label.appendChild(document.createTextNode(sig.name));
    checkboxContainer.appendChild(label);
  });

  // Select All Logic
  selectAllBox.addEventListener('change', (e) => {
    const isChecked = e.target.checked;

    // Update UI Boxes
    signalBoxes.forEach(box => box.checked = isChecked);

    // Update Plotly (Batch)
    const allIndices = signals.map(s => s.index);
    Plotly.restyle(plotDiv, { visible: isChecked }, allIndices);
  });
}

function downloadPlotly() {
  if (plotlyGraphDiv) {
    Plotly.downloadImage(plotlyGraphDiv, { format: 'png', filename: 'simulation_result' });
  }
}

function toggleSignalDropdown() {
  const content = document.getElementById('signal-dropdown-content');
  if (content) {
    content.classList.toggle('app-hidden');
  }
}

// Close dropdown when clicking outside
window.addEventListener('click', function (e) {
  const dropdown = document.querySelector('.custom-dropdown');
  const content = document.getElementById('signal-dropdown-content');
  const btn = document.getElementById('signal-dropdown-btn');

  if (dropdown && content && btn && !content.classList.contains('app-hidden')) {
    if (!dropdown.contains(e.target)) {
      content.classList.add('app-hidden');
    }
  }
});

// 4. Helper: Render Table
function renderTable(data) {
  const table = document.getElementById('data-table');
  const thead = table.querySelector('thead');
  const tbody = table.querySelector('tbody');

  // Clear existing
  thead.innerHTML = '';
  tbody.innerHTML = '';

  // Headers
  const trHead = document.createElement('tr');
  data.headers.forEach(header => {
    const th = document.createElement('th');
    th.textContent = header;
    trHead.appendChild(th);
  });
  thead.appendChild(trHead);

  // Rows (Limit to first 100 for performance if massive)
  const maxRows = 100; // Limit display to avoid DOM freeze
  data.rows.slice(0, maxRows).forEach(rowData => {
    const tr = document.createElement('tr');
    rowData.forEach(cell => {
      const td = document.createElement('td');
      // Format numbers to 4 decimal places if they are float numbers matches simple regex or typeof
      if (typeof cell === 'number') {
        td.textContent = cell.toPrecision(5);
      } else {
        td.textContent = cell;
      }
      tr.appendChild(td);
    });
    tbody.appendChild(tr);
  });

  // Notice if truncated
  if (data.rows.length > maxRows) {
    const tr = document.createElement('tr');
    const td = document.createElement('td');
    td.colSpan = data.headers.length;
    td.style.textAlign = 'center';
    td.style.fontStyle = 'italic';
    td.textContent = `... Showing first ${maxRows} of ${data.rows.length} rows. Download CSV for full data. ...`;
    tr.appendChild(td);
    tbody.appendChild(tr);
  }
}

// 5. Download Functions
function downloadPlot() {
  if (!globalPlotImage) return;
  const link = document.createElement('a');
  link.href = 'data:image/png;base64,' + globalPlotImage;
  link.download = 'mspice_simulation_plot.png';
  document.body.appendChild(link);
  link.click();
  document.body.removeChild(link);
}

function downloadCSV() {
  if (!globalCSVData) return;

  // Construct CSV content
  const headers = globalCSVData.headers.join(',');
  const rows = globalCSVData.rows.map(row => row.join(',')).join('\n');
  const csvContent = headers + '\n' + rows;

  // Create Blob
  const blob = new Blob([csvContent], { type: 'text/csv;charset=utf-8;' });
  const url = URL.createObjectURL(blob);

  // Trigger download
  const link = document.createElement('a');
  link.setAttribute('href', url);
  link.setAttribute('download', 'mspice_simulation_data.csv');
  link.style.visibility = 'hidden';
  document.body.appendChild(link);
  link.click();
  document.body.removeChild(link);
}