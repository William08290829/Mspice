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
      if (scrollDown <= sectionTop + sectionHeight){
        sectionsClass.classList.add('paw-active')
      } else {
        sectionsClass.classList.remove('paw-active')
      }
    } else {
      if (scrollDown > sectionTop && scrollDown <= sectionTop + sectionHeight){
        sectionsClass.classList.add('paw-active')
      } else {
        sectionsClass.classList.remove('paw-active')
      }
    }
  })
}

window.addEventListener('scroll', scrollActive)