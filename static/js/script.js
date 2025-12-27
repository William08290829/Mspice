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
    }, 1000);
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
async function runSimulation() {
  if (!netlistInput || !outputArea) return;

  const netlist = netlistInput.value;
  if (!netlist.trim()) {
    alert("Please enter or upload a netlist first!");
    return;
  }

  // Show loading state in output
  outputArea.innerHTML = '<div class="spinner" style="width:30px;height:30px;"></div><p>Simulating...</p>';

  try {
    const response = await fetch('/api/simulate', {
      method: 'POST',
      headers: {
        'Content-Type': 'application/json'
      },
      body: JSON.stringify({ netlist: netlist })
    });

    const data = await response.json();

    if (response.ok && data.image) {
      outputArea.innerHTML = `<img src="data:image/png;base64,${data.image}" alt="Simulation Result">`;
    } else if (data.error) {
      outputArea.innerHTML = `<p style="color:red;">Error: ${data.error}</p>`;
    } else {
      outputArea.innerHTML = `<p>${data.message || 'Unknown result'}</p>`;
    }

  } catch (err) {
    console.error(err);
    outputArea.innerHTML = `<p style="color:red;">Network Error: ${err.message}</p>`;
  }
}